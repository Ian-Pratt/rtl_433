/** @file
    Offline evaluator for the Honeywell CM921 soft FSK scanner profile.
*/

#include "soft_fsk.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint8_t *read_cu8(char const *path, size_t *sample_count)
{
    FILE *file = fopen(path, "rb");
    uint8_t *data;
    long length;

    if (!file) {
        fprintf(stderr, "%s: %s\n", path, strerror(errno));
        return NULL;
    }
    if (fseek(file, 0, SEEK_END) != 0 || (length = ftell(file)) < 2 ||
            fseek(file, 0, SEEK_SET) != 0) {
        fclose(file);
        return NULL;
    }
    *sample_count = (size_t)length / 2;
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
    static uint8_t const preamble[] = {0x55, 0xff, 0x00, 0x40};
    soft_fsk_result_t results[8];
    size_t sample_count;
    uint8_t *iq;
    unsigned output_bits;
    unsigned symbol_samples;

    if (argc != 6) {
        fprintf(stderr, "usage: %s FILE.cu8 TONE0_HZ TONE1_HZ SYMBOL_SAMPLES OUTPUT_BITS\n", argv[0]);
        return 2;
    }
    symbol_samples = strtoul(argv[4], NULL, 0);
    output_bits = strtoul(argv[5], NULL, 0);
    iq = read_cu8(argv[1], &sample_count);
    if (!iq || output_bits > SOFT_FSK_MAX_BYTES * 8) {
        free(iq);
        return 2;
    }

    int found = soft_fsk_scan_cu8_multi_framed(iq, sample_count, 1000000,
            strtof(argv[2], NULL), strtof(argv[3], NULL), symbol_samples,
            preamble, 26, output_bits, 26, 10, results,
            sizeof(results) / sizeof(results[0]));
    free(iq);
    if (found <= 0) {
        return 2;
    }

    for (int candidate = 0; candidate < found; ++candidate) {
        soft_fsk_result_t *result = &results[candidate];
        printf("score=%.3f start=%zu rate=%.3f hard=", result->score, result->start,
                result->symbol_period);
        for (unsigned i = 0; i < (result->bit_count + 7) / 8; ++i) {
            printf("%s%02x", i ? " " : "", result->bytes[i]);
        }
        putchar('\n');
    }
    return 0;
}
