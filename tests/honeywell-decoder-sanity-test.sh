#!/bin/sh
set -eu

rtl_433="$1"

# Fully framed demodulated bit rows: a normal 1f09 packet, a legitimate
# 2309/7eff binding packet, and a checksum-valid extended 3ef0 packet.
# Unusual but exactly framed, checksum-valid payloads must be retained for
# downstream analysis rather than rejected by an incomplete opcode allow-list.
normal='{306}55ff00599555954acb32acb34ccd5555334ad552ad2d2ad4b555552ad5354cd3354d5565555540'
binding='{306}55ff00599555954acb34cd554ccb52b4ab32d4b2ad2d2ad4b2acab54cd5555552ad33565555540'
extended_3ef0='{646}55ff00599555954acb332d2b4d4ab2ad4b52cd5554ab4acb32acab2ad2b2accb2ad4b2acb32ad332acd32ad532acad2ad2d2accd2ad4d2acb52ad352acd52ad554acab4ad2b4accb4ad4b554ad56555554'

normal_out=$($rtl_433 -R 0 -R 185 -F json -y "$normal" 2>/dev/null)
binding_out=$($rtl_433 -R 0 -R 185 -F json -y "$binding" 2>/dev/null)
extended_out=$($rtl_433 -R 0 -R 185 -F json -y "$extended_3ef0" 2>/dev/null)

printf '%s\n' "$normal_out" | grep -q '"Command" : "1f09"'
printf '%s\n' "$binding_out" | grep -q '"Command" : "2309"'
printf '%s\n' "$binding_out" | grep -q '"Payload" : "007eff"'
printf '%s\n' "$extended_out" | grep -q '"Command" : "3ef0"'
printf '%s\n' "$extended_out" | grep -q '"Payload" : "000102030405060708090a0b0c0d0e0f10111213"'
