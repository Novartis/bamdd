#!/usr/bin/env bash
# warm_xref_index.sh <src-dir> <qmd-basename> [<qmd-basename> ...]
#
# Pre-populates Quarto's project-scoped crossref index mapping file
# (<src-dir>/.quarto/xref/INDEX) with an entry for every <qmd, html/md/pdf>
# pair BEFORE the parallel `make -jN website` fan-out starts. All three
# output extensions are registered unconditionally (cheap - this never
# renders anything) because the shared INDEX key is (input, output
# basename), and the Makefile has %.html, %.md and %.pdf pattern rules that
# can all run in parallel against the same src/_quarto.yml project.
#
# Why: Quarto's crossrefIndexForOutputFile() (in quarto.js) only performs a
# non-atomic read-modify-write of the shared INDEX file the FIRST time a
# given (input, outputBaseFile) key is registered; once the key exists it
# only reads the file. Concurrent `quarto render` processes racing to
# register a *new* key concurrently corrupt INDEX (interleaved writes ->
# invalid JSON), which then breaks every subsequent render that tries to
# parse it. By registering all keys ourselves up front, in one fast serial
# pass (no rendering required), every parallel `quarto render` call only
# ever reads INDEX - which is safe under concurrency - so the expensive
# per-file R/Stan work still runs fully in parallel.
#
# This script is idempotent: it never touches or overwrites existing valid
# entries, it only adds missing ones, so re-running it (e.g. because a
# chapter file changed) is cheap and safe.
set -euo pipefail

if [ $# -lt 2 ]; then
  echo "Usage: $0 <src-dir> <qmd-basename> [<qmd-basename> ...]" >&2
  exit 1
fi

SRC_DIR=$1
shift

XREF_DIR="$SRC_DIR/.quarto/xref"
INDEX_FILE="$XREF_DIR/INDEX"
LOCK_FILE="$XREF_DIR/.warm.lock"

mkdir -p "$XREF_DIR"

# Zero-cost safety net: this whole script is meant to run once, serially,
# as a Make prerequisite, but flock guards against accidental concurrent
# invocations without needing any locking around the actual parallel
# `quarto render` calls themselves.
exec 200>"$LOCK_FILE"
flock 200

# Load the existing INDEX, tolerating a missing or already-corrupted file
# (this also makes the script self-healing if INDEX was left corrupt by a
# previous run).
json='{}'
if [ -s "$INDEX_FILE" ]; then
  if candidate=$(jq -c '.' "$INDEX_FILE" 2>/dev/null); then
    json=$candidate
  else
    echo "warm_xref_index.sh: $INDEX_FILE is not valid JSON, resetting it" >&2
  fi
fi

gen_id() {
  # 8 lowercase hex chars, matching the style of ids Quarto itself
  # generates (not required by Quarto - any unique string works - but kept
  # consistent for readability if you ever inspect INDEX by hand).
  local candidate
  while true; do
    candidate=$(od -An -tx1 -N4 /dev/urandom | tr -d ' \n')
    if ! grep -q "\"$candidate\"" <<<"$json"; then
      echo "$candidate"
      return
    fi
  done
}

changed=0
for qmd in "$@"; do
  base=${qmd%.qmd}
  for ext in html md pdf; do
    out="${base}.${ext}"
    existing=$(jq -r --arg q "$qmd" --arg o "$out" '.[$q][$o] // empty' <<<"$json")
    if [ -z "$existing" ]; then
      id=$(gen_id)
      json=$(jq -c --arg q "$qmd" --arg o "$out" --arg id "$id" \
        '.[$q] = ((.[$q] // {}) + {($o): $id})' <<<"$json")
      changed=1
    fi
  done
done

if [ "$changed" = "1" ] || [ ! -f "$INDEX_FILE" ]; then
  tmp=$(mktemp "$XREF_DIR/.INDEX.tmp.XXXXXX")
  jq '.' <<<"$json" > "$tmp"
  mv "$tmp" "$INDEX_FILE"
fi
