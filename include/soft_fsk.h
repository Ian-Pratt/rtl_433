/** @file
    Experimental narrow-channel soft FSK preamble scanner.
*/

#ifndef INCLUDE_SOFT_FSK_H_
#define INCLUDE_SOFT_FSK_H_

#include <stddef.h>
#include <stdint.h>

#define SOFT_FSK_MAX_BYTES 32
#define SOFT_FSK_STREAM_OVERLAP_SAMPLES 10000

typedef struct soft_fsk_result {
    size_t start;
    float score;
    unsigned symbol_samples;
    unsigned bit_count;
    uint8_t bytes[SOFT_FSK_MAX_BYTES];
} soft_fsk_result_t;

typedef struct soft_fsk_stream {
    uint8_t overlap[SOFT_FSK_STREAM_OVERLAP_SAMPLES * 2];
    size_t overlap_samples;
    uint64_t last_packet_start_plus_one;
} soft_fsk_stream_t;

/** Search CU8 IQ for an MSB-first preamble carried by two known FSK tones.

    The symbol rate and tones are deliberately explicit: this is a small
    protocol-specific fallback experiment, not another general demodulator.

    @return 1 when a candidate was found, 0 for insufficient input, or -1 on
    allocation/argument failure.
*/
int soft_fsk_scan_cu8(uint8_t const *iq, size_t sample_count, unsigned sample_rate,
        float tone0_hz, float tone1_hz, unsigned symbol_samples,
        uint8_t const *preamble, unsigned preamble_bits, unsigned output_bits,
        soft_fsk_result_t *result);

#endif /* INCLUDE_SOFT_FSK_H_ */
