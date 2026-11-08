#!/bin/sh
set -eu

rtl_433=${1:-"$(dirname "$0")/../src/rtl_433"}
temporary=$(mktemp -d)
trap 'rm -rf "$temporary"' EXIT HUP INT TERM
sample="$temporary/empty_868.3M_1000k.cu8"
stdout="$temporary/stdout"
stderr="$temporary/stderr"
: >"$sample"

"$rtl_433" -q -r "$sample" -f 868.3M -s 1000k -R 185 \
    -Y honeywellsoft -M statslog:3:60 -F json >"$stdout" 2>"$stderr"

test ! -s "$stdout"
line=$(grep '^rtl_433_stats [0-9][0-9]* {' "$stderr")
test -n "$line"
case "$line" in
    *'"soft_fsk":[{"profile":"honeywell"'*) ;;
    *) exit 1 ;;
esac
