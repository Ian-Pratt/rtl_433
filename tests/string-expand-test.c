#include "data.h"
#include "string_expand.h"

#include <stdio.h>
#include <string.h>

static char *no_sanitize(char *start, char *end)
{
    (void)start;
    return end;
}

static int check_expand(data_t *data, char const *format, char const *expected)
{
    char topic[256];
    expand_topic_string(topic, format, data, "receiver", no_sanitize);
    if (strcmp(topic, expected)) {
        fprintf(stderr, "FAIL: \"%s\" expanded to \"%s\", expected \"%s\"\n", format, topic, expected);
        return 1;
    }
    return 0;
}

int main(void)
{
    data_t *data = data_make(
            "model",   "", DATA_STRING, "Honeywell-CM921",
            "id",      "", DATA_INT,    158,
            "ids",     "", DATA_STRING, "0fefe9 0fefe9",
            "Command", "", DATA_STRING, "30c9",
            NULL);

    int failed = 0;
    failed += check_expand(data, "rtl_433[/model][/ids][/Command]", "rtl_433/Honeywell-CM921/0fefe9 0fefe9/30c9");
    failed += check_expand(data, "rtl_433/[Command]", "rtl_433/30c9");
    failed += check_expand(data, "rtl_433[/id][/missing]", "rtl_433/158");
    failed += check_expand(data, "rtl_433[/missing:fallback]", "rtl_433/fallback");
    failed += check_expand(data, "rtl_433[/hostname]", "rtl_433/receiver");

    data_free(data);
    return failed;
}
