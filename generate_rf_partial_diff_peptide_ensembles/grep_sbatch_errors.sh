#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <search_dir>" >&2
  exit 1
fi

search_dir="$1"

if [[ ! -d "$search_dir" ]]; then
  echo "Directory not found: $search_dir" >&2
  exit 1
fi

echo "Searching for 'Error' or 'CANCEL' in .out files under $search_dir..."
find "$search_dir" -type f -name "*.out" -print0 \
  | xargs -0 grep -lE "Error|CANCEL" \
  || echo "No matches found."