/** @file
    String formatting a-la MQTT topic

    Copyright (C) 2019 Christian Zuckschwerdt

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
*/

#include "string_expand.h"
#include "logger.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *append_topic(char *topic, data_t *data, expand_string_sanitizer sanitizer)
{
    if (data->type == DATA_STRING) {
        strcpy(topic, data->value.v_ptr); // NOLINT
        size_t len = strlen(topic);
        (*sanitizer)(topic, topic + len);
        topic += len;
    }
    else if (data->type == DATA_INT) {
        topic += sprintf(topic, "%d", data->value.v_int);
    }
    else {
        print_logf(LOG_ERROR, __func__, "Can't append data type %d to topic", data->type);
    }

    return topic;
}

static data_t *find_data_token(data_t *data, char const *start, char const *end)
{
    size_t len = end - start;
    for (data_t *d = data; d; d = d->next) {
        if (strlen(d->key) == len && !strncmp(d->key, start, len))
            return d;
    }
    return NULL;
}

char *expand_topic_string(char *topic, char const *format, data_t *data, char const *hostname, expand_string_sanitizer sanitizer)
{
    // consume entire format string
    while (format && *format) {
        data_t *data_token  = NULL;
        char const *string_token = NULL;
        int leading_slash   = 0;
        char const *t_start = NULL;
        char const *t_end   = NULL;
        char const *d_start = NULL;
        char const *d_end   = NULL;
        // copy until '['
        while (*format && *format != '[')
            *topic++ = *format++;
        // skip '['
        if (!*format)
            break;
        ++format;
        // read slash
        if (*format == '/') {
            leading_slash = *format;
            format++;
        }
        // read key until : or ]
        t_start = t_end = format;
        while (*format && *format != ':' && *format != ']' && *format != '[')
            t_end = ++format;
        // read default until ]
        if (*format == ':') {
            d_start = d_end = ++format;
            while (*format && *format != ']' && *format != '[')
                d_end = ++format;
        }
        // check for proper closing
        if (*format != ']') {
            print_log(LOG_FATAL, __func__, "unterminated token");
            exit(1);
        }
        ++format;

        // resolve token
        if (t_end - t_start == 8 && !strncmp(t_start, "hostname", 8))
            string_token = hostname;
        else
            data_token = find_data_token(data, t_start, t_end);

        // append token or default
        if (!data_token && !string_token && !d_start)
            continue;
        if (leading_slash)
            *topic++ = leading_slash;
        if (data_token)
            topic = append_topic(topic, data_token, sanitizer);
        else if (string_token)
            topic += sprintf(topic, "%s", string_token);
        else
            topic += sprintf(topic, "%.*s", (int)(d_end - d_start), d_start);
    }

    *topic = '\0';
    return topic;
}
