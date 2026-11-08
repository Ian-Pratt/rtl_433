/** @file
    Offline evaluator for the shared Fine Offset soft FSK scanner profile.
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
    uint8_t preamble[10];
    soft_fsk_result_t results[8];
    size_t sample_count;
    uint8_t *iq;
    unsigned symbol_samples;
    unsigned output_bytes;
    unsigned preamble_aa_bytes;

    if (argc != 7) {
        fprintf(stderr, "usage: %s FILE.cu8 TONE0_HZ TONE1_HZ SYMBOL_SAMPLES OUTPUT_BYTES PREAMBLE_AA_BYTES\n",
                argv[0]);
        return 2;
    }
    symbol_samples = strtoul(argv[4], NULL, 0);
    output_bytes = strtoul(argv[5], NULL, 0);
    preamble_aa_bytes = strtoul(argv[6], NULL, 0);
    iq = read_cu8(argv[1], &sample_count);
    if (!iq || output_bytes > SOFT_FSK_MAX_BYTES || preamble_aa_bytes > 8) {
        free(iq);
        return 2;
    }
    memset(preamble, 0xaa, preamble_aa_bytes);
    preamble[preamble_aa_bytes] = 0x2d;
    preamble[preamble_aa_bytes + 1] = 0xd4;
    int found = soft_fsk_scan_cu8_multi(iq, sample_count, 1000000,
                strtof(argv[2], NULL), strtof(argv[3], NULL), symbol_samples,
                preamble, (preamble_aa_bytes + 2) * 8, output_bytes * 8,
                results, sizeof(results) / sizeof(results[0]));
    if (found <= 0) {
        free(iq);
        return 2;
    }
    free(iq);

    for (int candidate = 0; candidate < found; ++candidate) {
        soft_fsk_result_t *result = &results[candidate];
        printf("score=%.3f start=%zu rate=%u hard=", result->score, result->start,
                result->symbol_samples);
        for (unsigned i = 0; i < result->bit_count / 8; ++i) {
            printf("%s%02x", i ? " " : "", result->bytes[i]);
        }
        putchar('\n');
    }
    return 0;
}
