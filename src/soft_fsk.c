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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct complex_d {
    double re;
    double im;
} complex_d_t;

static complex_d_t complex_add(complex_d_t a, complex_d_t b)
{
    complex_d_t out = {a.re + b.re, a.im + b.im};
    return out;
}

static complex_d_t complex_sub(complex_d_t a, complex_d_t b)
{
    complex_d_t out = {a.re - b.re, a.im - b.im};
    return out;
}

static complex_d_t complex_mul(complex_d_t a, complex_d_t b)
{
    complex_d_t out = {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
    return out;
}

static complex_d_t complex_scale(complex_d_t a, double scale)
{
    complex_d_t out = {a.re * scale, a.im * scale};
    return out;
}

static complex_d_t cu8_sample(uint8_t const *iq, size_t index)
{
    complex_d_t out = {(double)iq[index * 2] - 127.5, (double)iq[index * 2 + 1] - 127.5};
    return out;
}

static complex_d_t rectangular_dft(uint8_t const *iq, size_t start, unsigned width, double radians)
{
    complex_d_t sum = {0.0, 0.0};
    complex_d_t oscillator = {1.0, 0.0};
    complex_d_t step = {cos(-radians), sin(-radians)};

    for (unsigned i = 0; i < width; ++i) {
        sum = complex_add(sum, complex_mul(cu8_sample(iq, start + i), oscillator));
        oscillator = complex_mul(oscillator, step);
    }
    return sum;
}

static void advance_rectangular_dft(complex_d_t *sum, complex_d_t leaving,
        complex_d_t entering, complex_d_t inverse_step, complex_d_t end_step)
{
    *sum = complex_mul(inverse_step,
            complex_add(complex_sub(*sum, leaving), complex_mul(entering, end_step)));
}

static int tone_energy(uint8_t const *iq, size_t sample_count, unsigned sample_rate,
        float tone_hz, unsigned width, float *energy)
{
    size_t const valid = sample_count - width + 1;
    double const center = 2.0 * M_PI * tone_hz / sample_rate;
    double const side = 2.0 * M_PI / (width - 1);
    double radians[3] = {center, center - side, center + side};
    complex_d_t sum[3];
    complex_d_t inverse_step[3];
    complex_d_t end_step[3];

    for (unsigned band = 0; band < 3; ++band) {
        sum[band] = rectangular_dft(iq, 0, width, radians[band]);
        inverse_step[band].re = cos(radians[band]);
        inverse_step[band].im = sin(radians[band]);
        end_step[band].re = cos(-radians[band] * width);
        end_step[band].im = sin(-radians[band] * width);
    }

    for (size_t start = 0; start < valid; ++start) {
        complex_d_t hann = complex_sub(complex_scale(sum[0], 0.5),
                complex_scale(complex_add(sum[1], sum[2]), 0.25));
        energy[start] = (float)(hann.re * hann.re + hann.im * hann.im);

        if (start + 1 >= valid) {
            continue;
        }
        if ((start & 0xffff) == 0xffff) {
            for (unsigned band = 0; band < 3; ++band) {
                sum[band] = rectangular_dft(iq, start + 1, width, radians[band]);
            }
        }
        else {
            complex_d_t leaving = cu8_sample(iq, start);
            complex_d_t entering = cu8_sample(iq, start + width);
            for (unsigned band = 0; band < 3; ++band) {
                advance_rectangular_dft(&sum[band], leaving, entering,
                        inverse_step[band], end_step[band]);
            }
        }
    }
    return 0;
}

static unsigned get_bit(uint8_t const *bytes, unsigned bit)
{
    return (bytes[bit / 8] >> (7 - bit % 8)) & 1;
}

int soft_fsk_scan_cu8(uint8_t const *iq, size_t sample_count, unsigned sample_rate,
        float tone0_hz, float tone1_hz, unsigned symbol_samples,
        uint8_t const *preamble, unsigned preamble_bits, unsigned output_bits,
        soft_fsk_result_t *result)
{
    float *llr;
    size_t valid;
    size_t last_start;
    size_t best_start = 0;
    double template_mean = 0.0;
    double template_norm = 0.0;
    double best_score = -1.0;

    if (!iq || !preamble || !result || !sample_rate || symbol_samples < 2 ||
            !preamble_bits || output_bits > SOFT_FSK_MAX_BYTES * 8 ||
            sample_count < (size_t)output_bits * symbol_samples) {
        return -1;
    }

    valid = sample_count - symbol_samples + 1;
    llr = malloc(valid * sizeof(*llr));
    if (!llr) {
        return -1;
    }
    tone_energy(iq, sample_count, sample_rate, tone0_hz, symbol_samples, llr);
    for (size_t i = 0; i < valid; ++i) {
        llr[i] = logf(llr[i] + 1.0f);
    }

    {
        float *mark_energy = malloc(valid * sizeof(*mark_energy));
        if (!mark_energy) {
            free(llr);
            return -1;
        }
        tone_energy(iq, sample_count, sample_rate, tone1_hz, symbol_samples, mark_energy);
        for (size_t i = 0; i < valid; ++i) {
            llr[i] = logf(mark_energy[i] + 1.0f) - llr[i];
        }
        free(mark_energy);
    }

    for (unsigned bit = 0; bit < preamble_bits; ++bit) {
        template_mean += get_bit(preamble, bit) ? 1.0 : -1.0;
    }
    template_mean /= preamble_bits;
    for (unsigned bit = 0; bit < preamble_bits; ++bit) {
        double value = (get_bit(preamble, bit) ? 1.0 : -1.0) - template_mean;
        template_norm += value * value;
    }
    template_norm = sqrt(template_norm);

    last_start = valid - (size_t)(output_bits - 1) * symbol_samples;
    for (size_t start = 0; start < last_start; ++start) {
        double mean = 0.0;
        double norm = 0.0;
        double dot = 0.0;
        for (unsigned bit = 0; bit < preamble_bits; ++bit) {
            mean += llr[start + (size_t)bit * symbol_samples];
        }
        mean /= preamble_bits;
        for (unsigned bit = 0; bit < preamble_bits; ++bit) {
            double sample = llr[start + (size_t)bit * symbol_samples] - mean;
            double expected = (get_bit(preamble, bit) ? 1.0 : -1.0) - template_mean;
            dot += sample * expected;
            norm += sample * sample;
        }
        if (norm > 0.0) {
            double score = dot / (sqrt(norm) * template_norm);
            if (score > best_score) {
                best_score = score;
                best_start = start;
            }
        }
    }

    memset(result, 0, sizeof(*result));
    result->start = best_start;
    result->score = (float)best_score;
    result->symbol_samples = symbol_samples;
    result->bit_count = output_bits;
    for (unsigned bit = 0; bit < output_bits; ++bit) {
        if (llr[best_start + (size_t)bit * symbol_samples] > 0.0f) {
            result->bytes[bit / 8] |= 1 << (7 - bit % 8);
        }
    }
    free(llr);
    return 1;
}
