#!/usr/bin/env bash
# check_w_fast.sh  — prints files MISSING variable "w"
# Usage: ./check_w_fast.sh /path/to/dir ["*.nc"]

set -euo pipefail

dir="${1:-.}"
glob="${2:-*.nc}"

# treat empty globs as empty list
shopt -s nullglob

# header types we might see in ncdump -h
types='(byte|ubyte|char|string|short|ushort|int|uint|int64|uint64|float|double)'

missing_any=0
for f in "$dir"/$glob; do
  # look for lines like: "float w(..."
  if ! ncdump -h "$f" | grep -i "w(valid_time" > /dev/null  ; then
    echo "[MISSING] $f"
    missing_any=1
  fi
done

exit "$missing_any"
