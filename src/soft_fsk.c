/** @file
    Experimental narrow-channel soft FSK preamble scanner.

    A Hann-windowed detector is evaluated at every possible symbol start.  A
    Hann DFT is the sum of three rectangular DFTs, so each tone can be updated
    in constant time with a sliding recurrence rather than convolving every
    sample with a complete window.
*/

#include "soft_fsk.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#if defined(__ARM_NEON)
#include <arm_neon.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SOFT_FSK_PREFIX_BITS 24
#define SOFT_FSK_PREFIX_SCORE_SQ_MIN 0.16

typedef struct complex_f {
    float re;
    float im;
} complex_f_t;

static complex_f_t complex_add(complex_f_t a, complex_f_t b)
{
    complex_f_t out = {a.re + b.re, a.im + b.im};
    return out;
}

static complex_f_t complex_mul(complex_f_t a, complex_f_t b)
{
    complex_f_t out = {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
    return out;
}

static complex_f_t cu8_sample(uint8_t const *iq, size_t index)
{
    complex_f_t out = {(float)iq[index * 2] - 127.5f, (float)iq[index * 2 + 1] - 127.5f};
    return out;
}

static complex_f_t rectangular_dft(uint8_t const *iq, size_t start, unsigned width, float radians)
{
    complex_f_t sum = {0.0f, 0.0f};
    complex_f_t oscillator = {1.0f, 0.0f};
    complex_f_t step = {cosf(-radians), sinf(-radians)};

    for (unsigned i = 0; i < width; ++i) {
        sum = complex_add(sum, complex_mul(cu8_sample(iq, start + i), oscillator));
        oscillator = complex_mul(oscillator, step);
    }
    return sum;
}

static void advance_rectangular_dfts(float *sum_re, float *sum_im,
        float const *inverse_re, float const *inverse_im,
        float const *end_re, float const *end_im,
        complex_f_t leaving, complex_f_t entering)
{
#if defined(__ARM_NEON)
    float32x4_t leave_re = vdupq_n_f32(leaving.re);
    float32x4_t leave_im = vdupq_n_f32(leaving.im);

    for (unsigned band = 0; band < 8; band += 4) {
        float32x4_t inv_re = vld1q_f32(&inverse_re[band]);
        float32x4_t inv_im = vld1q_f32(&inverse_im[band]);
        float32x4_t finish_re = vld1q_f32(&end_re[band]);
        float32x4_t finish_im = vld1q_f32(&end_im[band]);
        float32x4_t enter_re = vfmsq_n_f32(vmulq_n_f32(finish_re, entering.re),
                finish_im, entering.im);
        float32x4_t enter_im = vfmaq_n_f32(vmulq_n_f32(finish_im, entering.re),
                finish_re, entering.im);
        float32x4_t shifted_re = vaddq_f32(vsubq_f32(vld1q_f32(&sum_re[band]), leave_re), enter_re);
        float32x4_t shifted_im = vaddq_f32(vsubq_f32(vld1q_f32(&sum_im[band]), leave_im), enter_im);
        float32x4_t next_re = vfmsq_f32(vmulq_f32(inv_re, shifted_re), inv_im, shifted_im);
        float32x4_t next_im = vfmaq_f32(vmulq_f32(inv_re, shifted_im), inv_im, shifted_re);
        vst1q_f32(&sum_re[band], next_re);
        vst1q_f32(&sum_im[band], next_im);
    }
#else
    for (unsigned band = 0; band < 8; ++band) {
        float enter_re = entering.re * end_re[band] - entering.im * end_im[band];
        float enter_im = entering.re * end_im[band] + entering.im * end_re[band];
        float shifted_re = sum_re[band] - leaving.re + enter_re;
        float shifted_im = sum_im[band] - leaving.im + enter_im;
        sum_re[band] = inverse_re[band] * shifted_re - inverse_im[band] * shifted_im;
        sum_im[band] = inverse_re[band] * shifted_im + inverse_im[band] * shifted_re;
    }
#endif
}

static int tone_pair_llr(uint8_t const *iq, size_t sample_count, unsigned sample_rate,
        float tone0_hz, float tone1_hz, unsigned width, float *llr)
{
    size_t const valid = sample_count - width + 1;
    float const side = 2.0f * (float)M_PI / (width - 1);
    float const center0 = 2.0f * (float)M_PI * tone0_hz / sample_rate;
    float const center1 = 2.0f * (float)M_PI * tone1_hz / sample_rate;
    /* Pad the six Hann resonators to two four-lane SIMD vectors. */
    float radians[8] = {
        center0, center0 - side, center0 + side,
        center1, center1 - side, center1 + side,
        center1, center1,
    };
    float sum_re[8];
    float sum_im[8];
    float inverse_re[8];
    float inverse_im[8];
    float end_re[8];
    float end_im[8];

    for (unsigned band = 0; band < 8; ++band) {
        complex_f_t sum = rectangular_dft(iq, 0, width, radians[band]);
        sum_re[band] = sum.re;
        sum_im[band] = sum.im;
        inverse_re[band] = cosf(radians[band]);
        inverse_im[band] = sinf(radians[band]);
        end_re[band] = cosf(-radians[band] * width);
        end_im[band] = sinf(-radians[band] * width);
    }

    for (size_t start = 0; start < valid; ++start) {
        float hann0_re = 0.5f * sum_re[0] - 0.25f * (sum_re[1] + sum_re[2]);
        float hann0_im = 0.5f * sum_im[0] - 0.25f * (sum_im[1] + sum_im[2]);
        float hann1_re = 0.5f * sum_re[3] - 0.25f * (sum_re[4] + sum_re[5]);
        float hann1_im = 0.5f * sum_im[3] - 0.25f * (sum_im[4] + sum_im[5]);
        float energy0 = hann0_re * hann0_re + hann0_im * hann0_im;
        float energy1 = hann1_re * hann1_re + hann1_im * hann1_im;
        llr[start] = (energy1 - energy0) / (energy1 + energy0 + 1.0f);

        if (start + 1 >= valid) {
            continue;
        }
        if ((start & 0xffff) == 0xffff) {
            for (unsigned band = 0; band < 8; ++band) {
                complex_f_t sum = rectangular_dft(iq, start + 1, width, radians[band]);
                sum_re[band] = sum.re;
                sum_im[band] = sum.im;
            }
        }
        else {
            complex_f_t leaving = cu8_sample(iq, start);
            complex_f_t entering = cu8_sample(iq, start + width);
            advance_rectangular_dfts(sum_re, sum_im, inverse_re, inverse_im,
                    end_re, end_im, leaving, entering);
        }
    }
    return 0;
}

static unsigned get_bit(uint8_t const *bytes, unsigned bit)
{
    return (bytes[bit / 8] >> (7 - bit % 8)) & 1;
}

static void score_candidate(float const *llr, size_t start, unsigned symbol_samples,
        double const *expected, unsigned preamble_bits, double template_norm_sq,
        soft_fsk_result_t *results, unsigned max_results, unsigned *result_count)
{
    double sum = 0.0;
    double sum_sq = 0.0;
    double dot = 0.0;

    for (unsigned bit = 0; bit < preamble_bits; ++bit) {
        double sample = llr[start + (size_t)bit * symbol_samples];
        sum += sample;
        sum_sq += sample * sample;
        dot += sample * expected[bit];
    }
    double norm_sq = sum_sq - sum * sum / preamble_bits;
    if (dot <= 0.0 || norm_sq <= 0.0) {
        return;
    }

    double score_sq = dot * dot / (norm_sq * template_norm_sq);
    unsigned match = *result_count;
    size_t separation = (size_t)preamble_bits * symbol_samples;
    for (unsigned i = 0; i < *result_count; ++i) {
        size_t distance = start > results[i].start ?
                start - results[i].start : results[i].start - start;
        if (distance < separation) {
            match = i;
            break;
        }
    }
    if (match < *result_count) {
        if (score_sq > results[match].score) {
            results[match].start = start;
            results[match].score = (float)score_sq;
        }
    }
    else if (*result_count < max_results) {
        results[*result_count].start = start;
        results[*result_count].score = (float)score_sq;
        ++*result_count;
    }
    else {
        unsigned weakest = 0;
        for (unsigned i = 1; i < *result_count; ++i) {
            if (results[i].score < results[weakest].score) {
                weakest = i;
            }
        }
        if (score_sq > results[weakest].score) {
            results[weakest].start = start;
            results[weakest].score = (float)score_sq;
        }
    }
}

int soft_fsk_scan_cu8_multi(uint8_t const *iq, size_t sample_count, unsigned sample_rate,
        float tone0_hz, float tone1_hz, unsigned symbol_samples,
        uint8_t const *preamble, unsigned preamble_bits, unsigned output_bits,
        soft_fsk_result_t *results, unsigned max_results)
{
    float *llr;
    size_t valid;
    size_t last_start;
    double template_mean = 0.0;
    double template_norm_sq = 0.0;
    double expected[SOFT_FSK_MAX_BYTES * 8];
    double prefix_mean = 0.0;
    double prefix_norm_sq = 0.0;
    double prefix_expected[SOFT_FSK_PREFIX_BITS];
#if defined(__ARM_NEON)
    float prefix_expected_f[SOFT_FSK_PREFIX_BITS];
#endif
    unsigned prefix_bits;
    unsigned result_count = 0;

    if (!iq || !preamble || !results || !max_results || !sample_rate || symbol_samples < 2 ||
            !preamble_bits || output_bits > SOFT_FSK_MAX_BYTES * 8 ||
            sample_count < (size_t)output_bits * symbol_samples) {
        return -1;
    }
    memset(results, 0, max_results * sizeof(*results));

    valid = sample_count - symbol_samples + 1;
    llr = malloc(valid * sizeof(*llr));
    if (!llr) {
        return -1;
    }
    tone_pair_llr(iq, sample_count, sample_rate, tone0_hz, tone1_hz, symbol_samples, llr);

    for (unsigned bit = 0; bit < preamble_bits; ++bit) {
        template_mean += get_bit(preamble, bit) ? 1.0 : -1.0;
    }
    template_mean /= preamble_bits;
    for (unsigned bit = 0; bit < preamble_bits; ++bit) {
        expected[bit] = (get_bit(preamble, bit) ? 1.0 : -1.0) - template_mean;
        template_norm_sq += expected[bit] * expected[bit];
    }

    prefix_bits = preamble_bits < SOFT_FSK_PREFIX_BITS ? preamble_bits : SOFT_FSK_PREFIX_BITS;
    for (unsigned bit = 0; bit < prefix_bits; ++bit) {
        prefix_mean += get_bit(preamble, bit) ? 1.0 : -1.0;
    }
    prefix_mean /= prefix_bits;
    for (unsigned bit = 0; bit < prefix_bits; ++bit) {
        prefix_expected[bit] = (get_bit(preamble, bit) ? 1.0 : -1.0) - prefix_mean;
        prefix_norm_sq += prefix_expected[bit] * prefix_expected[bit];
#if defined(__ARM_NEON)
        prefix_expected_f[bit] = (float)prefix_expected[bit];
#endif
    }

    last_start = valid - (size_t)(output_bits - 1) * symbol_samples;
#if defined(__ARM_NEON)
    size_t start = 0;
    float prefix_limit = (float)(SOFT_FSK_PREFIX_SCORE_SQ_MIN * prefix_norm_sq);
    for (; start + 4 <= last_start; start += 4) {
        float32x4_t prefix_sum = vdupq_n_f32(0.0f);
        float32x4_t prefix_sum_sq = vdupq_n_f32(0.0f);
        float32x4_t prefix_dot = vdupq_n_f32(0.0f);
        for (unsigned bit = 0; bit < prefix_bits; ++bit) {
            float32x4_t sample = vld1q_f32(&llr[start + (size_t)bit * symbol_samples]);
            prefix_sum = vaddq_f32(prefix_sum, sample);
            prefix_sum_sq = vfmaq_f32(prefix_sum_sq, sample, sample);
            prefix_dot = vfmaq_n_f32(prefix_dot, sample, prefix_expected_f[bit]);
        }
        float32x4_t prefix_sample_norm_sq = vfmsq_n_f32(prefix_sum_sq,
                vmulq_f32(prefix_sum, prefix_sum), 1.0f / prefix_bits);
        float norms[4];
        float dots[4];
        vst1q_f32(norms, prefix_sample_norm_sq);
        vst1q_f32(dots, prefix_dot);
        for (unsigned lane = 0; lane < 4; ++lane) {
            if (dots[lane] > 0.0f && norms[lane] > 0.0f &&
                    dots[lane] * dots[lane] > prefix_limit * norms[lane]) {
                score_candidate(llr, start + lane, symbol_samples, expected,
                        preamble_bits, template_norm_sq, results, max_results, &result_count);
            }
        }
    }
    for (; start < last_start; ++start) {
#else
    for (size_t start = 0; start < last_start; ++start) {
#endif
        double prefix_sum = 0.0;
        double prefix_sum_sq = 0.0;
        double prefix_dot = 0.0;
        for (unsigned bit = 0; bit < prefix_bits; ++bit) {
            double sample = llr[start + (size_t)bit * symbol_samples];
            prefix_sum += sample;
            prefix_sum_sq += sample * sample;
            prefix_dot += sample * prefix_expected[bit];
        }
        double prefix_sample_norm_sq = prefix_sum_sq - prefix_sum * prefix_sum / prefix_bits;
        if (prefix_dot <= 0.0 || prefix_sample_norm_sq <= 0.0 ||
                prefix_dot * prefix_dot <= SOFT_FSK_PREFIX_SCORE_SQ_MIN *
                        prefix_sample_norm_sq * prefix_norm_sq) {
            continue;
        }

        score_candidate(llr, start, symbol_samples, expected, preamble_bits,
                template_norm_sq, results, max_results, &result_count);
    }

    for (unsigned i = 1; i < result_count; ++i) {
        soft_fsk_result_t candidate = results[i];
        unsigned j = i;
        while (j && results[j - 1].start > candidate.start) {
            results[j] = results[j - 1];
            --j;
        }
        results[j] = candidate;
    }
    for (unsigned i = 0; i < result_count; ++i) {
        size_t packet_start = results[i].start;
        float score_sq = results[i].score;
        memset(results[i].bytes, 0, sizeof(results[i].bytes));
        results[i].score = sqrtf(score_sq);
        results[i].symbol_samples = symbol_samples;
        results[i].bit_count = output_bits;
        for (unsigned bit = 0; bit < output_bits; ++bit) {
            if (llr[packet_start + (size_t)bit * symbol_samples] > 0.0f) {
                results[i].bytes[bit / 8] |= 1 << (7 - bit % 8);
            }
        }
    }
    free(llr);
    return result_count;
}

int soft_fsk_scan_cu8(uint8_t const *iq, size_t sample_count, unsigned sample_rate,
        float tone0_hz, float tone1_hz, unsigned symbol_samples,
        uint8_t const *preamble, unsigned preamble_bits, unsigned output_bits,
        soft_fsk_result_t *result)
{
    return soft_fsk_scan_cu8_multi(iq, sample_count, sample_rate, tone0_hz, tone1_hz,
            symbol_samples, preamble, preamble_bits, output_bits, result, 1);
}
