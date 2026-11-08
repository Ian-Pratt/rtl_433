/** @file
    Experimental narrow-channel soft FSK preamble scanner.
*/

#ifndef INCLUDE_SOFT_FSK_H_
#define INCLUDE_SOFT_FSK_H_

#include <stddef.h>
#include <stdint.h>

#define SOFT_FSK_MAX_BYTES 256
#define SOFT_FSK_STREAM_OVERLAP_SAMPLES 65536
#define SOFT_FSK_MAX_DEVICES 4

typedef struct soft_fsk_stats {
    unsigned scan_calls;
    unsigned candidates;
    unsigned below_score;
    unsigned unknown_family;
    unsigned filtered;
    unsigned incomplete;
    unsigned normal_suppressed;
    unsigned duplicate_suppressed;
    unsigned decode_attempts;
    unsigned decode_ok;
    unsigned fail_mic;
    unsigned fail_other;
} soft_fsk_stats_t;

typedef struct soft_fsk_result {
    size_t start;
    float score;
    float symbol_period;
    unsigned symbol_samples;
    unsigned bit_count;
    uint8_t bytes[SOFT_FSK_MAX_BYTES];
} soft_fsk_result_t;

typedef struct soft_fsk_stream {
    uint8_t overlap[SOFT_FSK_STREAM_OVERLAP_SAMPLES * 2];
    size_t overlap_samples;
    uint64_t last_packet_start_plus_one;
    uint64_t normal_block_begin[SOFT_FSK_MAX_DEVICES];
    uint64_t normal_block_end[SOFT_FSK_MAX_DEVICES];
    soft_fsk_stats_t stats;
} soft_fsk_stream_t;

/** Search CU8 IQ for an MSB-first preamble carried by two known FSK tones.

    The symbol rate and tones are deliberately explicit so protocol profiles
    can share the narrow-channel detector without a costly frequency search.

    @return 1 when a candidate was found, 0 for insufficient input, or -1 on
    allocation/argument failure.
*/
int soft_fsk_scan_cu8(uint8_t const *iq, size_t sample_count, unsigned sample_rate,
        float tone0_hz, float tone1_hz, unsigned symbol_samples,
        uint8_t const *preamble, unsigned preamble_bits, unsigned output_bits,
        soft_fsk_result_t *result);

/** As soft_fsk_scan_cu8(), but retain the strongest non-overlapping candidates.

    Results are ordered by sample position. This lets protocol decoders reject a
    damaged burst and try a repeated transmission without another DSP pass.
*/
int soft_fsk_scan_cu8_multi(uint8_t const *iq, size_t sample_count, unsigned sample_rate,
        float tone0_hz, float tone1_hz, unsigned symbol_samples,
        uint8_t const *preamble, unsigned preamble_bits, unsigned output_bits,
        soft_fsk_result_t *results, unsigned max_results);

/** As soft_fsk_scan_cu8_multi(), with bounded fractional clock recovery using
    asynchronous start/stop framing. frame_offset is the first start bit and
    frame_bits includes the zero start and one stop bit. */
int soft_fsk_scan_cu8_multi_framed(uint8_t const *iq, size_t sample_count,
        unsigned sample_rate, float tone0_hz, float tone1_hz, unsigned symbol_samples,
        uint8_t const *preamble, unsigned preamble_bits, unsigned output_bits,
        unsigned frame_offset, unsigned frame_bits,
        soft_fsk_result_t *results, unsigned max_results);

#endif /* INCLUDE_SOFT_FSK_H_ */
