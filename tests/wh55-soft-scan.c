/** @file
    Offline evaluator for the experimental WH55 soft FSK scanner.
*/

#include "soft_fsk.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint8_t crc8(uint8_t const *message, unsigned length)
{
    uint8_t remainder = 0;

    for (unsigned byte = 0; byte < length; ++byte) {
        remainder ^= message[byte];
        for (unsigned bit = 0; bit < 8; ++bit) {
            remainder = remainder & 0x80 ? (remainder << 1) ^ 0x31 : remainder << 1;
        }
    }
    return remainder;
}

static uint8_t *read_cu8(char const *path, size_t start_sample, size_t requested_samples,
        size_t *sample_count)
{
    FILE *file = fopen(path, "rb");
    uint8_t *data;
    long length;

    if (!file) {
        fprintf(stderr, "%s: %s\n", path, strerror(errno));
        return NULL;
    }
    if (fseek(file, 0, SEEK_END) != 0 || (length = ftell(file)) < 2) {
        fclose(file);
        return NULL;
    }
    if (start_sample > (size_t)length / 2) {
        fclose(file);
        return NULL;
    }
    *sample_count = (size_t)length / 2 - start_sample;
    if (requested_samples && requested_samples < *sample_count) {
        *sample_count = requested_samples;
    }
    if (fseek(file, (long)(start_sample * 2), SEEK_SET) != 0) {
        fclose(file);
        return NULL;
    }
    data = malloc(*sample_count * 2);
    if (!data || fread(data, 2, *sample_count, file) != *sample_count) {
        free(data);
        fclose(file);
        return NULL;
    }
    fclose(file);
    return data;
}

int main(int argc, char **argv)
{
    static uint8_t const preamble[] = {0xaa, 0xaa, 0xaa, 0x2d, 0xd4, 0x55};
    soft_fsk_result_t result;
    size_t sample_count;
    uint8_t *iq;
    float tone0;
    float tone1;
    size_t start_sample = 0;
    size_t requested_samples = 0;
    int valid_crc;

    if (argc < 4 || argc > 6) {
        fprintf(stderr, "usage: %s FILE.cu8 TONE0_HZ TONE1_HZ [START_SAMPLE [SAMPLE_COUNT]]\n",
                argv[0]);
        return 2;
    }
    tone0 = strtof(argv[2], NULL);
    tone1 = strtof(argv[3], NULL);
    if (argc > 4) {
        start_sample = strtoull(argv[4], NULL, 0);
    }
    if (argc > 5) {
        requested_samples = strtoull(argv[5], NULL, 0);
    }
    iq = read_cu8(argv[1], start_sample, requested_samples, &sample_count);
    if (!iq) {
        return 2;
    }
    if (soft_fsk_scan_cu8(iq, sample_count, 1000000, tone0, tone1, 58,
                preamble, 48, 17 * 8, &result) <= 0) {
        free(iq);
        return 2;
    }
    free(iq);

    printf("score=%.3f start=%zu rate=%u hard=", result.score,
            start_sample + result.start,
            result.symbol_samples);
    for (unsigned i = 0; i < result.bit_count / 8; ++i) {
        printf("%s%02x", i ? " " : "", result.bytes[i]);
    }
    valid_crc = result.bytes[5] == 0x55 &&
            crc8(&result.bytes[5], 8) == result.bytes[13];
    printf(" crc=%s\n", valid_crc ? "OK" : "BAD");
    return valid_crc ? 0 : 1;
}
