#!/bin/sh
set -e
MODE="$1"
LOG="tests/ci/${MODE}.log"
case "$MODE" in
  tap)
    sed -i 's/{ 1, 3, 0}/{ 1, 4, 0}/' globals.c
    ;;
  crc)
    perl -0pi -e 's/overall \^= __builtin_parity\(word & ~\(1u<<\(32-1\)\)\);/overall ^= __builtin_parity(word & ~(1u<<(32-1))) ^ 0x01;/' tests/helpers.c
    ;;
  *)
    echo "Usage: $0 tap|crc" >&2
    exit 1
    ;;
esac
if make check > "$LOG" 2>&1; then
  cat "$LOG"
  echo "Unexpected pass" >&2
  exit 1
else
  cat "$LOG"
  exit 0
fi
