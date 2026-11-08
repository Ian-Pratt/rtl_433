/** @file
    Generic RF data receiver and decoder for ISM band devices using RTL-SDR and SoapySDR.

    Copyright (C) 2026 Christian W. Zuckschwerdt <zany@triq.net>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "r_flow.h"
#include "rtl_433.h"
#include "r_private.h"
#include "r_device.h"
#include "r_api.h"
#include "bitbuffer.h"
#include "baseband.h"
#include "pulse_analyzer.h"
#include "pulse_detect.h"
#include "pulse_detect_fsk.h"
#include "pulse_slicer.h"
#include "data.h"
#include "raw_output.h"
#include "r_util.h"
#include "am_analyze.h"
#include "logger.h"
#include "fatal.h"
#include "soft_fsk.h"

#define FINEOFFSET_SOFT_TONE0_HZ 868320000.0f
#define FINEOFFSET_SOFT_TONE1_HZ 868375000.0f
#define FINEOFFSET_SOFT_SCORE_MIN 0.90f
#define FINEOFFSET_SOFT_MAX_CANDIDATES 8
#define HONEYWELL_SOFT_TONE0_HZ 868199000.0f
#define HONEYWELL_SOFT_TONE1_HZ 868338000.0f
#define HONEYWELL_SOFT_SCORE_MIN 0.70f
#define HONEYWELL_SOFT_MAX_CANDIDATES 8
#define HONEYWELL_SOFT_OUTPUT_BITS 1536

static r_device *find_device(struct dm_state *demod, char const *name)
{
    for (void **iter = demod->r_devs.elems; iter && *iter; ++iter) {
        r_device *device = *iter;
        if (!strcmp(device->name, name)) {
            return device;
        }
    }
    return NULL;
}

static float cu8_mean_power(uint8_t const *iq, size_t begin, size_t end)
{
    double power = 0.0;

    if (end <= begin) {
        return 1.0f;
    }
    for (size_t sample = begin; sample < end; ++sample) {
        float i = iq[sample * 2] - 127.5f;
        float q = iq[sample * 2 + 1] - 127.5f;
        power += i * i + q * q;
    }
    return (float)(power / (end - begin));
}

static float cu8_noise_floor(uint8_t const *iq, size_t sample_count,
        size_t packet_begin, size_t packet_end)
{
    float floor = INFINITY;
    size_t before = packet_begin > 20000 ? packet_begin - 20000 : 0;
    size_t after = packet_end + 20000 < sample_count ? packet_end + 20000 : sample_count;

    for (size_t begin = before; begin + 500 <= packet_begin; begin += 500) {
        float power = cu8_mean_power(iq, begin, begin + 500);
        floor = power < floor ? power : floor;
    }
    for (size_t begin = packet_end; begin + 500 <= after; begin += 500) {
        float power = cu8_mean_power(iq, begin, begin + 500);
        floor = power < floor ? power : floor;
    }
    return isfinite(floor) ? floor : 1.0f;
}

static int run_fineoffset_soft(r_cfg_t *cfg, r_device *wh55_device,
        r_device *wn34_device, r_device *wh31_device, r_device *wh51_device, uint8_t const *iq,
        size_t sample_count, int normal_wh55_decoded, int normal_wn34_decoded,
        int normal_wh31_decoded, int normal_wh51_decoded)
{
    static uint8_t const wh55_preamble[] = {0xaa, 0xaa, 0xaa, 0x2d, 0xd4, 0x55};
    static uint8_t const fineoffset_preamble[] = {0xaa, 0xaa, 0x2d, 0xd4};
    struct dm_state *demod = cfg->demod;
    soft_fsk_stream_t *stream = &demod->fineoffset_soft;
    soft_fsk_result_t results[FINEOFFSET_SOFT_MAX_CANDIDATES];
    size_t old_overlap = stream->overlap_samples;
    size_t combined_samples = old_overlap + sample_count;
    uint8_t *combined;
    unsigned symbol_samples;
    uint8_t const *preamble;
    unsigned preamble_bits;
    unsigned output_bits;
    unsigned family_index;
    float tone0;
    float tone1;
    int events = 0;

    if ((!cfg->wh55_soft_id && !cfg->fineoffset_soft) ||
            (!wh55_device && !wn34_device && !wh31_device && !wh51_device) ||
            demod->sample_size != 2 || !sample_count ||
            demod->samp_rate < 800000 || demod->samp_rate > 1200000) {
        return 0;
    }

    tone0 = FINEOFFSET_SOFT_TONE0_HZ - demod->center_frequency;
    tone1 = FINEOFFSET_SOFT_TONE1_HZ - demod->center_frequency;
    if (fabsf(tone0) >= demod->samp_rate / 2.0f || fabsf(tone1) >= demod->samp_rate / 2.0f) {
        return 0;
    }
    symbol_samples = (unsigned)lroundf(demod->samp_rate * 58.0e-6f);
    if (cfg->fineoffset_soft) {
        preamble = fineoffset_preamble;
        preamble_bits = sizeof(fineoffset_preamble) * 8;
        output_bits = 23 * 8;
        family_index = 4;
    }
    else {
        preamble = wh55_preamble;
        preamble_bits = sizeof(wh55_preamble) * 8;
        output_bits = 17 * 8;
        family_index = 5;
    }
    combined = malloc(combined_samples * 2);
    if (!combined) {
        return 0;
    }
    memcpy(combined, stream->overlap, old_overlap * 2);
    memcpy(combined + old_overlap * 2, iq, sample_count * 2);

    stream->stats.scan_calls++;
    int found = soft_fsk_scan_cu8_multi(combined, combined_samples, demod->samp_rate,
            tone0, tone1, symbol_samples, preamble, preamble_bits, output_bits,
            results, FINEOFFSET_SOFT_MAX_CANDIDATES);

    size_t keep = combined_samples < SOFT_FSK_STREAM_OVERLAP_SAMPLES ?
            combined_samples : SOFT_FSK_STREAM_OVERLAP_SAMPLES;
    memcpy(stream->overlap, combined + (combined_samples - keep) * 2, keep * 2);
    stream->overlap_samples = keep;

    if (found <= 0) {
        free(combined);
        return 0;
    }
    stream->stats.candidates += found;

    pulse_data_t saved_fsk = demod->fsk_pulse_data;
    unsigned saved_start_ago = demod->pulse_data.start_ago;
    uint64_t combined_pos = demod->input_pos - old_overlap;
    int normal_decoded[4] = {
        normal_wh55_decoded, normal_wn34_decoded, normal_wh31_decoded, normal_wh51_decoded,
    };
    for (unsigned i = 0; i < 4; ++i) {
        if (normal_decoded[i]) {
            stream->normal_block_begin[i] = combined_pos;
            stream->normal_block_end[i] = demod->input_pos + sample_count;
        }
    }
    for (int candidate = 0; candidate < found; ++candidate) {
        soft_fsk_result_t *result = &results[candidate];
        uint8_t family = result->bytes[family_index];
        r_device *device = NULL;
        unsigned device_index = 0;
        if (result->score < FINEOFFSET_SOFT_SCORE_MIN) {
            stream->stats.below_score++;
            continue;
        }
        if (family == 0x55 && wh55_device) {
            unsigned sensor_id = (result->bytes[family_index + 2] << 8) |
                    result->bytes[family_index + 3];
            if (cfg->wh55_soft_id && cfg->wh55_soft_id != UINT32_MAX &&
                            sensor_id != cfg->wh55_soft_id) {
                stream->stats.filtered++;
                continue;
            }
            device = wh55_device;
        }
        else if ((family == 0x34 || family == 0x38) && wn34_device) {
            device = wn34_device;
            device_index = 1;
        }
        else if ((family == 0x20 || family == 0x30 || family == 0x37 ||
                    family == 0x40 || family == 0x52 || family == 0x68) && wh31_device) {
            device = wh31_device;
            device_index = 2;
        }
        else if (family == 0x51 && wh51_device) {
            device = wh51_device;
            device_index = 3;
        }
        if (!device) {
            stream->stats.unknown_family++;
            continue;
        }
        if (normal_decoded[device_index]) {
            stream->stats.normal_suppressed++;
            continue;
        }

        size_t packet_samples = result->bit_count * result->symbol_samples;
        size_t settle_samples = (size_t)lroundf(demod->samp_rate * 3.0e-3f);
        if (result->start + packet_samples + settle_samples > combined_samples) {
            stream->stats.incomplete++;
            continue;
        }

        uint64_t packet_pos = combined_pos + result->start;
        if (stream->normal_block_end[device_index] &&
                packet_pos >= stream->normal_block_begin[device_index] &&
                packet_pos <= stream->normal_block_end[device_index]) {
            stream->stats.normal_suppressed++;
            continue;
        }

        if (stream->last_packet_start_plus_one &&
                packet_pos <= stream->last_packet_start_plus_one - 1 + symbol_samples) {
            stream->stats.duplicate_suppressed++;
            continue;
        }
        bitbuffer_t bits = {0};
        for (unsigned bit = 0; bit < result->bit_count; ++bit) {
            bitbuffer_add_bit(&bits, result->bytes[bit / 8] >> (7 - bit % 8) & 1);
        }

        size_t packet_end = result->start + packet_samples < combined_samples ?
                result->start + packet_samples : combined_samples;
        float signal_power = cu8_mean_power(combined, result->start, packet_end);
        float noise_power = cu8_noise_floor(combined, combined_samples, result->start, packet_end);

        memset(&demod->fsk_pulse_data, 0, sizeof(demod->fsk_pulse_data));
        demod->fsk_pulse_data.fsk_f1_est = 1;
        demod->fsk_pulse_data.fsk_f2_est = 1;
        demod->fsk_pulse_data.freq1_hz = FINEOFFSET_SOFT_TONE0_HZ;
        demod->fsk_pulse_data.freq2_hz = FINEOFFSET_SOFT_TONE1_HZ;
        demod->fsk_pulse_data.centerfreq_hz = demod->center_frequency;
        demod->fsk_pulse_data.rssi_db = 10.0f * log10f(signal_power) - 42.1442f;
        demod->fsk_pulse_data.noise_db = 10.0f * log10f(noise_power) - 42.1442f;
        demod->fsk_pulse_data.snr_db = 10.0f * log10f(signal_power / noise_power);
        demod->pulse_data.start_ago = packet_pos < demod->input_pos + sample_count ?
                demod->input_pos + sample_count - packet_pos : 0;

        stream->stats.decode_attempts++;
        int ret = device->decode_fn(device, &bits);
        device->decode_events += 1;
        if (ret > 0) {
            stream->last_packet_start_plus_one = packet_pos + 1;
            device->decode_ok += 1;
            device->decode_messages += ret;
            stream->stats.decode_ok += ret;
            events += ret;
        }
        else if (ret >= DECODE_FAIL_SANITY) {
            device->decode_fails[-ret] += 1;
            if (ret == DECODE_FAIL_MIC) {
                stream->stats.fail_mic++;
            }
            else {
                stream->stats.fail_other++;
            }
        }
    }

    demod->pulse_data.start_ago = saved_start_ago;
    demod->fsk_pulse_data = saved_fsk;
    free(combined);
    return events;
}

static int run_honeywell_soft(r_cfg_t *cfg, r_device *device, uint8_t const *iq,
        size_t sample_count, int normal_decoded)
{
    static uint8_t const preamble[] = {0x55, 0xff, 0x00, 0x40};
    struct dm_state *demod = cfg->demod;
    soft_fsk_stream_t *stream = &demod->honeywell_soft;
    soft_fsk_result_t results[HONEYWELL_SOFT_MAX_CANDIDATES];
    size_t old_overlap = stream->overlap_samples;
    size_t combined_samples = old_overlap + sample_count;
    uint8_t *combined;
    float tone0;
    float tone1;
    unsigned symbol_samples;
    int events = 0;

    if (!cfg->honeywell_soft || !device || demod->sample_size != 2 || !sample_count ||
            demod->samp_rate < 800000 || demod->samp_rate > 1200000) {
        return 0;
    }

    tone0 = HONEYWELL_SOFT_TONE0_HZ - demod->center_frequency;
    tone1 = HONEYWELL_SOFT_TONE1_HZ - demod->center_frequency;
    if (fabsf(tone0) >= demod->samp_rate / 2.0f || fabsf(tone1) >= demod->samp_rate / 2.0f) {
        return 0;
    }
    symbol_samples = (unsigned)lroundf(demod->samp_rate * 26.0e-6f);
    combined = malloc(combined_samples * 2);
    if (!combined) {
        return 0;
    }
    memcpy(combined, stream->overlap, old_overlap * 2);
    memcpy(combined + old_overlap * 2, iq, sample_count * 2);

    stream->stats.scan_calls++;
    int found = soft_fsk_scan_cu8_multi_framed(combined, combined_samples, demod->samp_rate,
            tone0, tone1, symbol_samples, preamble, 26, HONEYWELL_SOFT_OUTPUT_BITS,
            26, 10,
            results, HONEYWELL_SOFT_MAX_CANDIDATES);

    size_t keep = combined_samples < SOFT_FSK_STREAM_OVERLAP_SAMPLES ?
            combined_samples : SOFT_FSK_STREAM_OVERLAP_SAMPLES;
    memcpy(stream->overlap, combined + (combined_samples - keep) * 2, keep * 2);
    stream->overlap_samples = keep;

    if (found <= 0) {
        free(combined);
        return 0;
    }
    stream->stats.candidates += found;

    pulse_data_t saved_fsk = demod->fsk_pulse_data;
    unsigned saved_start_ago = demod->pulse_data.start_ago;
    uint64_t combined_pos = demod->input_pos - old_overlap;
    for (int candidate = 0; candidate < found; ++candidate) {
        soft_fsk_result_t *result = &results[candidate];
        if (result->score < HONEYWELL_SOFT_SCORE_MIN) {
            stream->stats.below_score++;
            continue;
        }
        if (normal_decoded) {
            stream->stats.normal_suppressed++;
            continue;
        }

        size_t packet_samples = (size_t)lroundf(result->bit_count * result->symbol_period);
        if (result->start + packet_samples > combined_samples) {
            stream->stats.incomplete++;
            continue;
        }
        uint64_t packet_pos = combined_pos + result->start;
        if (stream->last_packet_start_plus_one &&
                packet_pos <= stream->last_packet_start_plus_one - 1 + symbol_samples) {
            stream->stats.duplicate_suppressed++;
            continue;
        }

        bitbuffer_t bits = {0};
        for (unsigned bit = 0; bit < result->bit_count; ++bit) {
            bitbuffer_add_bit(&bits, result->bytes[bit / 8] >> (7 - bit % 8) & 1);
        }

        size_t packet_end = result->start + packet_samples;
        float signal_power = cu8_mean_power(combined, result->start, packet_end);
        float noise_power = cu8_noise_floor(combined, combined_samples, result->start, packet_end);
        memset(&demod->fsk_pulse_data, 0, sizeof(demod->fsk_pulse_data));
        demod->fsk_pulse_data.fsk_f1_est = 1;
        demod->fsk_pulse_data.fsk_f2_est = 1;
        demod->fsk_pulse_data.freq1_hz = HONEYWELL_SOFT_TONE0_HZ;
        demod->fsk_pulse_data.freq2_hz = HONEYWELL_SOFT_TONE1_HZ;
        demod->fsk_pulse_data.centerfreq_hz = demod->center_frequency;
        demod->fsk_pulse_data.rssi_db = 10.0f * log10f(signal_power) - 42.1442f;
        demod->fsk_pulse_data.noise_db = 10.0f * log10f(noise_power) - 42.1442f;
        demod->fsk_pulse_data.snr_db = 10.0f * log10f(signal_power / noise_power);
        demod->pulse_data.start_ago = packet_pos < demod->input_pos + sample_count ?
                demod->input_pos + sample_count - packet_pos : 0;

        stream->stats.decode_attempts++;
        int ret = device->decode_fn(device, &bits);
        device->decode_events += 1;
        if (ret > 0) {
            stream->last_packet_start_plus_one = packet_pos + 1;
            device->decode_ok += 1;
            device->decode_messages += ret;
            stream->stats.decode_ok += ret;
            events += ret;
        }
        else if (ret >= DECODE_FAIL_SANITY) {
            device->decode_fails[-ret] += 1;
            if (ret == DECODE_FAIL_MIC) {
                stream->stats.fail_mic++;
            }
            else {
                stream->stats.fail_other++;
            }
        }
    }

    demod->pulse_data.start_ago = saved_start_ago;
    demod->fsk_pulse_data = saved_fsk;
    free(combined);
    return events;
}

static void calc_rssi_snr(struct dm_state const *demod, pulse_data_t *pulse_data)
{
    float ook_high_estimate      = pulse_data->ook_high_estimate > 0 ? pulse_data->ook_high_estimate : 1;
    float ook_low_estimate       = pulse_data->ook_low_estimate > 0 ? pulse_data->ook_low_estimate : 1;
    int const OOK_MAX_HIGH_LEVEL = DB_TO_AMP(0); // Maximum estimate for high level (-0 dB)
    float ook_max_estimate       = ook_high_estimate < OOK_MAX_HIGH_LEVEL ? ook_high_estimate : OOK_MAX_HIGH_LEVEL;
    float asnr                   = ook_max_estimate / ook_low_estimate;
    float foffs1                 = (float)pulse_data->fsk_f1_est / INT16_MAX * demod->samp_rate / 2.0f;
    float foffs2                 = (float)pulse_data->fsk_f2_est / INT16_MAX * demod->samp_rate / 2.0f;
    pulse_data->freq1_hz         = (foffs1 + demod->center_frequency);
    pulse_data->freq2_hz         = (foffs2 + demod->center_frequency);
    pulse_data->centerfreq_hz    = demod->center_frequency;
    pulse_data->depth_bits       = demod->sample_size * 4;
    // NOTE: for (CU8) amplitude is 10x (because it's squares)
    if (demod->sample_size == 2 && !demod->use_mag_est) {                    // amplitude (CU8)
        pulse_data->range_db = 42.1442f;                                     // 10*log10f(16384.0f) == 20*log10f(128.0f)
        pulse_data->rssi_db  = 10.0f * log10f(ook_high_estimate) - 42.1442f; // 10*log10f(16384.0f)
        pulse_data->noise_db = 10.0f * log10f(ook_low_estimate) - 42.1442f;  // 10*log10f(16384.0f)
        pulse_data->snr_db   = 10.0f * log10f(asnr);
    }
    else {                               // magnitude (CU8, CS16)
        pulse_data->range_db = 84.2884f; // 20*log10f(16384.0f)
        // lowest (scaled x128) reading at  8 bit is -20*log10(128) = -42.1442 (eff. -36 dB)
        // lowest (scaled div2) reading at 12 bit is -20*log10(1024) = -60.2060 (eff. -54 dB)
        // lowest (scaled div2) reading at 16 bit is -20*log10(16384) = -84.2884 (eff. -78 dB)
        pulse_data->rssi_db  = 20.0f * log10f(ook_high_estimate) - 84.2884f; // 20*log10f(16384.0f)
        pulse_data->noise_db = 20.0f * log10f(ook_low_estimate) - 84.2884f;  // 20*log10f(16384.0f)
        pulse_data->snr_db   = 20.0f * log10f(asnr);
    }
}

/**
Flush the SDR IQ data frame processing, e.g. on the end of a input file.

@return Count of successful decoding events
*/
int flush_sdr_flow(r_cfg_t *cfg)
{
    return push_sdr_flow(cfg, NULL, 0);
}

/**
Reset the SDR IQ data frame processing, e.g. on a new input file.
*/
void reset_sdr_flow(r_cfg_t *cfg)
{
    struct dm_state *demod = cfg->demod;

    get_time_now(&demod->now);

    demod->frame_start_ago   = 0;
    demod->frame_end_ago     = 0;
    demod->frame_event_count = 0;
    demod->frame_quality     = 0;

    demod->min_level_auto = 0.0f;
    demod->noise_level    = 0.0f;
    memset(&demod->fineoffset_soft, 0, sizeof(demod->fineoffset_soft));
    memset(&demod->honeywell_soft, 0, sizeof(demod->honeywell_soft));

    baseband_low_pass_filter_reset(&demod->lowpass_filter_state);
    baseband_demod_FM_reset(&demod->demod_FM_state);

    pulse_detect_reset(demod->pulse_detect);
}

/**
Push an IQ data frame to the SDR IQ data frame processing.

@return Count of successful decoding events
*/
int push_sdr_flow(r_cfg_t *cfg, unsigned char *iq_buf, uint32_t len)
{
    //fprintf(stderr, "push_sdr_flow... %u\n", len);
    struct dm_state *demod = cfg->demod;
    char time_str[LOCAL_TIME_BUFLEN];

    if (!demod) {
        // might happen when the demod closed and we get a last data frame
        return 0; // ignore the data
    }

    unsigned long n_samples = len / demod->sample_size;
    int soft_enabled = cfg->wh55_soft_id || cfg->fineoffset_soft || cfg->honeywell_soft;
    r_device *wh55_device = soft_enabled ? find_device(demod,
            "Fine Offset / Ecowitt WH55 water leak sensor") : NULL;
    r_device *wn34_device = cfg->fineoffset_soft ? find_device(demod,
            "Fine Offset Electronics WN34S/L/D, WN38 and Froggit DP150/D35 temperature sensor") : NULL;
    r_device *wh31_device = cfg->fineoffset_soft ? find_device(demod,
            "Ambient Weather WH31E Thermo-Hygrometer Sensor, EcoWitt WH40/WN20 rain gauge, WS68 weather station") : NULL;
    r_device *wh51_device = cfg->fineoffset_soft ? find_device(demod,
            "Fine Offset Electronics/Ecowitt WH51, WN31, SwitchDoc Labs SM23 Soil Moisture Sensor") : NULL;
    r_device *honeywell_device = cfg->honeywell_soft ? find_device(demod,
            "Honeywell CM921 Wireless Programmable Room Thermostat") : NULL;
    unsigned wh55_decode_ok_before = wh55_device ? wh55_device->decode_ok : 0;
    unsigned wn34_decode_ok_before = wn34_device ? wn34_device->decode_ok : 0;
    unsigned wh31_decode_ok_before = wh31_device ? wh31_device->decode_ok : 0;
    unsigned wh51_decode_ok_before = wh51_device ? wh51_device->decode_ok : 0;
    unsigned honeywell_decode_ok_before = honeywell_device ? honeywell_device->decode_ok : 0;
    if (n_samples * demod->sample_size != len) {
        print_log(LOG_WARNING, __func__, "Sample buffer length not aligned to sample size!");
    }

    int process_frame = 1;

    // Process new frame data if available
    if (len) {

    // Feed data to all raw outputs (e.g. rtl_tcp)
    // do this here and not in sdr_handler so realtime replay can use rtl_tcp output
    for (void **iter = demod->raw_handler->elems; iter && *iter; ++iter) {
        raw_output_t *output = *iter;
        raw_output_frame(output, iq_buf, len);
    }

    // save last frame time to see if a new second started
    time_t last_frame_sec = demod->now.tv_sec;
    get_time_now(&demod->now);

    // age the frame position if there is one
    if (demod->frame_start_ago) {
        demod->frame_start_ago += n_samples;
    }
    if (demod->frame_end_ago) {
        demod->frame_end_ago += n_samples;
    }

    // Feed data to the sample grabber if available
    if (demod->samp_grab) {
        samp_grab_push(demod->samp_grab, iq_buf, len);
    }

    // AM demodulation
    float avg_db;
    if (demod->sample_size == 2) { // CU8
        if (demod->use_mag_est) {
            //magnitude_true_cu8(iq_buf, demod->buf.temp, n_samples);
            avg_db = magnitude_est_cu8(iq_buf, demod->buf.temp, n_samples);
        }
        else { // amp est
            avg_db = envelope_detect(iq_buf, demod->buf.temp, n_samples);
        }
    } else { // CS16
        //magnitude_true_cs16((int16_t *)iq_buf, demod->buf.temp, n_samples);
        avg_db = magnitude_est_cs16((int16_t *)iq_buf, demod->buf.temp, n_samples);
    }

    // Squelch silent frames
    //fprintf(stderr, "noise level: %.1f dB current: %.1f dB min level: %.1f dB\n", demod->noise_level, avg_db, demod->min_level_auto);
    if (demod->min_level_auto == 0.0f) {
        demod->min_level_auto = demod->min_level;
    }
    if (demod->noise_level == 0.0f) {
        demod->noise_level = demod->min_level_auto - 3.0f;
    }
    int noise_only = avg_db < demod->noise_level + 3.0f; // or demod->min_level_auto?
    // always process frames if loader, dumper, or analyzers are in use, otherwise skip silent frames
    process_frame = demod->squelch_offset <= 0 || !noise_only || demod->load_info.format || demod->analyze_pulses || demod->dumper.len || demod->samp_grab;
    demod->total_frames_count += 1;
    if (noise_only) {
        demod->total_frames_squelch += 1;
        demod->noise_level = (demod->noise_level * 7 + avg_db) / 8; // fast fall over 8 frames
        // If auto_level and noise level well below min_level and significant change in noise level
        if (demod->auto_level > 0 && demod->noise_level < demod->min_level - 3.0f
                && fabsf(demod->min_level_auto - demod->noise_level - 3.0f) > 1.0f) {
            demod->min_level_auto = demod->noise_level + 3.0f;
            print_logf(LOG_WARNING, "Auto Level", "Estimated noise level is %.1f dB, adjusting minimum detection level to %.1f dB",
                    demod->noise_level, demod->min_level_auto);
            pulse_detect_set_levels(demod->pulse_detect, demod->use_mag_est, demod->level_limit, demod->min_level_auto, demod->min_snr, demod->detect_verbosity);
        }
    } else {
        demod->noise_level = (demod->noise_level * 31 + avg_db) / 32; // slow rise over 32 frames
    }
    // Report noise every report_noise seconds, but only for the first frame that second
    if (demod->report_noise && last_frame_sec != demod->now.tv_sec && demod->now.tv_sec % demod->report_noise == 0) {
        print_logf(LOG_WARNING, "Auto Level", "Current %s level %.1f dB, estimated noise %.1f dB",
                noise_only ? "noise" : "signal", avg_db, demod->noise_level);
    }

    // Run AM filters
    if (process_frame) {
        baseband_low_pass_filter(&demod->lowpass_filter_state, demod->buf.temp, demod->am_buf, n_samples);
    }

    // FM demodulation
    if (demod->enable_FM_demod && process_frame) {
        // Use manual FM low pass value otherwise 0.2 for minmax or 0.1 for classic
        float fm_low_pass = demod->fm_low_pass != 0.0f ? demod->fm_low_pass : demod->fsk_pulse_detect_mode ? 0.2f : 0.1f;
        if (demod->sample_size == 2) { // CU8
            baseband_demod_FM(&demod->demod_FM_state, iq_buf, demod->buf.fm, n_samples, demod->samp_rate, fm_low_pass);
        } else { // CS16
            baseband_demod_FM_cs16(&demod->demod_FM_state, (int16_t *)iq_buf, demod->buf.fm, n_samples, demod->samp_rate, fm_low_pass);
        }
    }

    // Handle special input formats
    if (demod->load_info.format == S16_AM) { // The IQ buffer is really AM demodulated data
        if (len > sizeof(demod->am_buf)) {
            FATAL("Buffer too small");
        }
        memcpy(demod->am_buf, iq_buf, len);
    } else if (demod->load_info.format == S16_FM) { // The IQ buffer is really FM demodulated data
        // we would need AM for the envelope too
        if (len > sizeof(demod->buf.fm)) {
            FATAL("Buffer too small");
        }
        memcpy(demod->buf.fm, iq_buf, len);
    }

    }

    // Run a pulse discriminator and pass packages to all configured slicers
    int d_events = 0; // Sensor events successfully detected
    if (demod->r_devs.len || demod->analyze_pulses || demod->dumper.len || demod->samp_grab) {
        // Detect a package and loop through demodulators with pulse data
        int package_type = PULSE_DATA_OOK;  // Just to get us started
        // Initialize all U8 logic buffers
        for (void **iter = demod->dumper.elems; iter && *iter; ++iter) {
            file_info_t const *dumper = *iter;
            if (dumper->format == U8_LOGIC) {
                memset(demod->u8_buf, 0, n_samples);
                break;
            }
        }
        while (package_type && process_frame) {
            int p_events = 0; // Sensor events successfully detected per package
            package_type = pulse_detect_package(demod->pulse_detect, demod->am_buf, demod->buf.fm, n_samples,
                    demod->samp_rate, demod->input_pos, &demod->pulse_data, &demod->fsk_pulse_data,
                    &demod->fsk_pulse_data_classic, demod->fsk_pulse_detect_mode);
            if (package_type) {
                // new package: set a first frame start if we are not tracking one already
                if (!demod->frame_start_ago) {
                    demod->frame_start_ago = demod->pulse_data.start_ago;
                }
                // always update the last frame end
                demod->frame_end_ago = demod->pulse_data.end_ago;
            }
            if (package_type == PULSE_DATA_OOK) {
                calc_rssi_snr(demod, &demod->pulse_data);
                if (demod->analyze_pulses) {
                    fprintf(stderr, "Detected OOK package\t%s\n", time_pos_str(cfg, demod->pulse_data.start_ago, time_str));
                }

                p_events += run_ook_demods(&demod->r_devs, &demod->pulse_data);
                demod->total_frames_ook += 1;
                demod->total_frames_events += p_events > 0;
                demod->frames_ook += 1;
                demod->frames_events += p_events > 0;

                // Dump pulse data for this complete package
                for (void **iter = demod->dumper.elems; iter && *iter; ++iter) {
                    file_info_t const *dumper = *iter;
                    if (dumper->format == VCD_LOGIC) {
                        pulse_data_print_vcd(dumper->file, &demod->pulse_data, '\'');
                    }
                    if (dumper->format == U8_LOGIC) {
                        pulse_data_dump_raw(demod->u8_buf, n_samples, demod->input_pos, &demod->pulse_data, 0x02);
                    }
                    if (dumper->format == PULSE_OOK) {
                        pulse_data_dump(dumper->file, &demod->pulse_data);
                    }
                }

                if (demod->verbosity >= LOG_TRACE) {
                    pulse_data_print(&demod->pulse_data);
                }
                if (demod->raw_mode == 1 || (demod->raw_mode == 2 && p_events == 0) || (demod->raw_mode == 3 && p_events > 0)) {
                    data_t *data = pulse_data_print_data(&demod->pulse_data);
                    event_occurred_handler(cfg, data);
                }
                if (demod->analyze_pulses && (demod->grab_mode <= 1 || (demod->grab_mode == 2 && p_events == 0) || (demod->grab_mode == 3 && p_events > 0))) {
                    r_device device = {.log_fn = log_device_handler, .output_ctx = cfg};
                    pulse_analyzer(&demod->pulse_data, package_type, &device);
                }
                if (demod->grab_mode == 4 && p_events == 0) {
                    r_device device = {.log_fn = log_device_handler, .output_ctx = cfg};
                    int p_quality   = pulse_analyzer_check(&demod->pulse_data, package_type, &device);
                    demod->frame_quality = p_quality > demod->frame_quality ? p_quality : demod->frame_quality;
                }
            }
            else if (package_type == PULSE_DATA_FSK) {
                calc_rssi_snr(demod, &demod->fsk_pulse_data);
                if (demod->analyze_pulses) {
                    fprintf(stderr, "Detected FSK package\t%s\n", time_pos_str(cfg, demod->fsk_pulse_data.start_ago, time_str));
                }

                p_events += run_fsk_demods(&demod->r_devs, &demod->fsk_pulse_data);
                if (p_events == 0 && demod->fsk_pulse_detect_mode == FSK_PULSE_DETECT_HYBRID) {
                    pulse_data_t fsk_pulse_data = demod->fsk_pulse_data;
                    demod->fsk_pulse_data = demod->fsk_pulse_data_classic;
                    demod->fsk_pulse_data_classic = fsk_pulse_data;
                    calc_rssi_snr(demod, &demod->fsk_pulse_data);
                    p_events += run_fsk_demods(&demod->r_devs, &demod->fsk_pulse_data);
                }
                demod->total_frames_fsk += 1;
                demod->total_frames_events += p_events > 0;
                demod->frames_fsk += 1;
                demod->frames_events += p_events > 0;

                // Dump pulse data for this complete package
                for (void **iter = demod->dumper.elems; iter && *iter; ++iter) {
                    file_info_t const *dumper = *iter;
                    if (dumper->format == VCD_LOGIC) {
                        pulse_data_print_vcd(dumper->file, &demod->fsk_pulse_data, '"');
                    }
                    if (dumper->format == U8_LOGIC) {
                        pulse_data_dump_raw(demod->u8_buf, n_samples, demod->input_pos, &demod->fsk_pulse_data, 0x04);
                    }
                    if (dumper->format == PULSE_OOK) {
                        pulse_data_dump(dumper->file, &demod->fsk_pulse_data);
                    }
                }

                if (demod->verbosity >= LOG_TRACE) {
                    pulse_data_print(&demod->fsk_pulse_data);
                }
                if (demod->raw_mode == 1 || (demod->raw_mode == 2 && p_events == 0) || (demod->raw_mode == 3 && p_events > 0)) {
                    data_t *data = pulse_data_print_data(&demod->fsk_pulse_data);
                    event_occurred_handler(cfg, data);
                }
                if (demod->analyze_pulses && (demod->grab_mode <= 1 || (demod->grab_mode == 2 && p_events == 0) || (demod->grab_mode == 3 && p_events > 0))) {
                    r_device device = {.log_fn = log_device_handler, .output_ctx = cfg};
                    pulse_analyzer(&demod->fsk_pulse_data, package_type, &device);
                }
                if (demod->grab_mode == 4 && p_events == 0) {
                    r_device device = {.log_fn = log_device_handler, .output_ctx = cfg};
                    int p_quality   = pulse_analyzer_check(&demod->fsk_pulse_data, package_type, &device);
                    demod->frame_quality = p_quality > demod->frame_quality ? p_quality : demod->frame_quality;
                }
            } // if (package_type == ...
            d_events += p_events;
        } // while (package_type)...

        // add event counter to the frames currently tracked
        demod->frame_event_count += d_events;

        // end frame tracking if older than a whole buffer
        if (demod->frame_start_ago && demod->frame_end_ago > n_samples) {
            if (demod->samp_grab) {
                if (demod->grab_mode == 1
                        || (demod->grab_mode == 2 && demod->frame_event_count == 0)
                        || (demod->grab_mode == 3 && demod->frame_event_count > 0)
                        || (demod->grab_mode == 4 && demod->frame_event_count == 0 && demod->frame_quality > 0)) {
                    unsigned frame_pad = n_samples / 8; // this could also be a fixed value, e.g. 10000 samples
                    unsigned start_padded = demod->frame_start_ago + frame_pad;
                    unsigned end_padded = demod->frame_end_ago - frame_pad;
                    unsigned len_padded = start_padded - end_padded;
                    samp_grab_write(demod->samp_grab, len_padded, end_padded);
                }
            }
            demod->frame_start_ago   = 0;
            demod->frame_event_count = 0;
            demod->frame_quality     = 0;
        }

        // Dump partial pulse data, might overlap with the last complete package
        for (void **iter = demod->dumper.elems; iter && *iter; ++iter) {
            file_info_t const *dumper = *iter;
            if (dumper->format == U8_LOGIC) {
                pulse_data_dump_raw(demod->u8_buf, n_samples, demod->input_pos, &demod->pulse_data, 0x02);
                pulse_data_dump_raw(demod->u8_buf, n_samples, demod->input_pos, &demod->fsk_pulse_data, 0x04);
                break;
            }
        }
    }

    // End processing if no new frame data
    if (!len) {
        return d_events;
    }

    // Run the AM analyzer (deprecated)
    if (demod->am_analyze) {
        am_analyze(demod->am_analyze, demod->am_buf, n_samples, demod->verbosity >= LOG_INFO, NULL);
    }

    // Save data to all dumpers (expect logic dumpers)
    for (void **iter = demod->dumper.elems; iter && *iter; ++iter) {
        file_info_t const *dumper = *iter;
        if (!dumper->file
                || dumper->format == VCD_LOGIC
                || dumper->format == PULSE_OOK) {
            continue;
        }
        uint8_t *out_buf = iq_buf;  // Default is to dump IQ samples
        unsigned long out_len = n_samples * demod->sample_size;

        if (dumper->format == CU8_IQ) {
            if (demod->sample_size == 4) {
                for (unsigned long n = 0; n < n_samples * 2; ++n)
                    ((uint8_t *)demod->buf.temp)[n] = (((int16_t *)iq_buf)[n] / 256) + 128; // scale Q0.15 to Q0.7
                out_buf = (uint8_t *)demod->buf.temp;
                out_len = n_samples * 2 * sizeof(uint8_t);
            }
        }
        else if (dumper->format == CS16_IQ) {
            if (demod->sample_size == 2) {
                for (unsigned long n = 0; n < n_samples * 2; ++n)
                    ((int16_t *)demod->buf.temp)[n] = (iq_buf[n] * 256) - 32768; // scale Q0.7 to Q0.15
                out_buf = (uint8_t *)demod->buf.temp; // this buffer is too small if out_block_size is large
                out_len = n_samples * 2 * sizeof(int16_t);
            }
        }
        else if (dumper->format == CS8_IQ) {
            if (demod->sample_size == 2) {
                for (unsigned long n = 0; n < n_samples * 2; ++n)
                    ((int8_t *)demod->buf.temp)[n] = (iq_buf[n] - 128);
            }
            else if (demod->sample_size == 4) {
                for (unsigned long n = 0; n < n_samples * 2; ++n)
                    ((int8_t *)demod->buf.temp)[n] = ((int16_t *)iq_buf)[n] >> 8;
            }
            out_buf = (uint8_t *)demod->buf.temp;
            out_len = n_samples * 2 * sizeof(int8_t);
        }
        else if (dumper->format == CF32_IQ) {
            if (demod->sample_size == 2) {
                for (unsigned long n = 0; n < n_samples * 2; ++n)
                    ((float *)demod->buf.temp)[n] = (iq_buf[n] - 128) / 128.0f;
            }
            else if (demod->sample_size == 4) {
                for (unsigned long n = 0; n < n_samples * 2; ++n)
                    ((float *)demod->buf.temp)[n] = ((int16_t *)iq_buf)[n] / 32768.0f;
            }
            out_buf = (uint8_t *)demod->buf.temp; // this buffer is too small if out_block_size is large
            out_len = n_samples * 2 * sizeof(float);
        }
        else if (dumper->format == S16_AM) {
            out_buf = (uint8_t *)demod->am_buf;
            out_len = n_samples * sizeof(int16_t);
        }
        else if (dumper->format == S16_FM) {
            out_buf = (uint8_t *)demod->buf.fm;
            out_len = n_samples * sizeof(int16_t);
        }
        else if (dumper->format == F32_AM) {
            for (unsigned long n = 0; n < n_samples; ++n)
                demod->f32_buf[n] = demod->am_buf[n] * (1.0f / 0x8000); // scale from Q0.15
            out_buf = (uint8_t *)demod->f32_buf;
            out_len = n_samples * sizeof(float);
        }
        else if (dumper->format == F32_FM) {
            for (unsigned long n = 0; n < n_samples; ++n)
                demod->f32_buf[n] = demod->buf.fm[n] * (1.0f / 0x8000); // scale from Q0.15
            out_buf = (uint8_t *)demod->f32_buf;
            out_len = n_samples * sizeof(float);
        }
        else if (dumper->format == F32_I) {
            if (demod->sample_size == 2) {
                for (unsigned long n = 0; n < n_samples; ++n)
                    demod->f32_buf[n] = (iq_buf[n * 2] - 128) * (1.0f / 0x80); // scale from Q0.7
            }
            else {
                for (unsigned long n = 0; n < n_samples; ++n)
                    demod->f32_buf[n] = ((int16_t *)iq_buf)[n * 2] * (1.0f / 0x8000); // scale from Q0.15
            }
            out_buf = (uint8_t *)demod->f32_buf;
            out_len = n_samples * sizeof(float);
        }
        else if (dumper->format == F32_Q) {
            if (demod->sample_size == 2) {
                for (unsigned long n = 0; n < n_samples; ++n)
                    demod->f32_buf[n] = (iq_buf[n * 2 + 1] - 128) * (1.0f / 0x80); // scale from Q0.7
            }
            else {
                for (unsigned long n = 0; n < n_samples; ++n)
                    demod->f32_buf[n] = ((int16_t *)iq_buf)[n * 2 + 1] * (1.0f / 0x8000); // scale from Q0.15
            }
            out_buf = (uint8_t *)demod->f32_buf;
            out_len = n_samples * sizeof(float);
        }
        else if (dumper->format == U8_LOGIC) { // state data
            out_buf = demod->u8_buf;
            out_len = n_samples;
        }

        if (fwrite(out_buf, 1, out_len, dumper->file) != out_len) {
            print_log(LOG_ERROR, __func__, "Short write, samples lost, exiting!");
            d_events = -1; // report error, handler should exit
        }
    }

    if (len && (wh55_device || wn34_device || wh31_device || wh51_device)) {
        d_events += run_fineoffset_soft(cfg, wh55_device, wn34_device, wh31_device, wh51_device,
                iq_buf, n_samples,
                wh55_device && wh55_device->decode_ok > wh55_decode_ok_before,
                wn34_device && wn34_device->decode_ok > wn34_decode_ok_before,
                wh31_device && wh31_device->decode_ok > wh31_decode_ok_before,
                wh51_device && wh51_device->decode_ok > wh51_decode_ok_before);
    }
    if (len && honeywell_device) {
        d_events += run_honeywell_soft(cfg, honeywell_device, iq_buf, n_samples,
                honeywell_device->decode_ok > honeywell_decode_ok_before);
    }

    demod->input_pos += n_samples;

    return d_events;
}
