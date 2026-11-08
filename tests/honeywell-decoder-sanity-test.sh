#!/bin/sh
set -eu

rtl_433="$1"

# Fully framed demodulated bit rows: a normal 1f09 packet, a legitimate
# 2309/7eff binding packet, and a checksum-valid 3ef0 with an impossible
# 20-byte payload. The last row used to be emitted as a valid event.
normal='{306}55ff00599555954acb32acb34ccd5555334ad552ad2d2ad4b555552ad5354cd3354d5565555540'
binding='{306}55ff00599555954acb34cd554ccb52b4ab32d4b2ad2d2ad4b2acab54cd5555552ad33565555540'
bad_shape='{646}55ff00599555954acb332d2b4d4ab2ad4b52cd5554ab4acb32acab2ad2b2accb2ad4b2acb32ad332acd32ad532acad2ad2d2accd2ad4d2acb52ad352acd52ad554acab4ad2b4accb4ad4b554ad56555554'

normal_out=$($rtl_433 -R 0 -R 185 -F json -y "$normal" 2>/dev/null)
binding_out=$($rtl_433 -R 0 -R 185 -F json -y "$binding" 2>/dev/null)
bad_out=$($rtl_433 -R 0 -R 185 -F json -y "$bad_shape" 2>/dev/null || true)

printf '%s\n' "$normal_out" | grep -q '"Command" : "1f09"'
printf '%s\n' "$binding_out" | grep -q '"Command" : "2309"'
printf '%s\n' "$binding_out" | grep -q '"Payload" : "007eff"'
test -z "$bad_out"
