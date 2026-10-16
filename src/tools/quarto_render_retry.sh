#!/usr/bin/env bash
# quarto_render_retry.sh <project-dir> [--] <command...>
#
# Defense-in-depth wrapper around a `quarto render` invocation. The primary
# fix for the parallel crossref INDEX race is warm_xref_index.sh, which
# pre-registers every (qmd, html/md/pdf) key before the parallel `make -jN`
# fan-out starts, so this path should rarely trigger. It exists to make the
# build self-healing against cases the pre-warm step doesn't cover (e.g. an
# ad-hoc render with a custom --output name, or a newly added .qmd file
# that slipped in between the warm-up pass and the render).
#
# If the wrapped command fails with output matching Quarto's known
# crossref-index JSON-parse race signature, this script RESTORES
# <project-dir>/.quarto/xref/INDEX by re-running warm_xref_index.sh over
# every *.qmd file found in <project-dir> (rather than just deleting
# INDEX): a bare deletion would wipe out every other in-flight render's
# already-registered key at once, which - under -jN - can hand several
# still-running `quarto render` processes a now-missing key at the same
# moment and reproduce the very write-write race this whole fix exists to
# prevent. Re-warming instead restores the same fully-populated, read-only
# state the build started with. Falls back to a plain `rm -f` only if the
# warm-up script itself can't be found/run. Retries a bounded number of
# times; any other failure is passed straight through unmodified (exit
# code and output preserved), so real render errors are never masked.
set -uo pipefail

if [ $# -lt 2 ]; then
  echo "Usage: $0 <project-dir> [--] <command...>" >&2
  exit 1
fi

PROJECT_DIR=$1
shift
if [ "${1:-}" = "--" ]; then
  shift
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WARM_SCRIPT="$SCRIPT_DIR/warm_xref_index.sh"

XREF_DIR="$PROJECT_DIR/.quarto/xref"
INDEX_FILE="$XREF_DIR/INDEX"

CROSSREF_ERROR_PATTERN='crossrefIndexForOutputFile|Unexpected non-whitespace character after JSON'
MAX_ATTEMPTS=3

repair_index() {
  mkdir -p "$XREF_DIR"
  if [ -x "$WARM_SCRIPT" ]; then
    local qmds=()
    while IFS= read -r -d '' f; do
      qmds+=("$(basename "$f")")
    done < <(find "$PROJECT_DIR" -maxdepth 1 -name '*.qmd' -print0 2>/dev/null)
    if [ "${#qmds[@]}" -gt 0 ] && "$WARM_SCRIPT" "$PROJECT_DIR" "${qmds[@]}"; then
      return 0
    fi
  fi
  # Fallback: warm-up script unavailable/failed - just drop the corrupted
  # file so Quarto recreates a fresh, valid one on the next attempt.
  rm -f "$INDEX_FILE"
}

attempt=1
while true; do
  out=$("$@" 2>&1)
  status=$?

  if [ $status -eq 0 ]; then
    printf '%s\n' "$out"
    exit 0
  fi

  if [ $attempt -ge $MAX_ATTEMPTS ] || ! grep -qE "$CROSSREF_ERROR_PATTERN" <<<"$out"; then
    printf '%s\n' "$out" >&2
    exit $status
  fi

  echo "quarto_render_retry.sh: detected crossref INDEX race (attempt $attempt/$MAX_ATTEMPTS), restoring $INDEX_FILE and retrying..." >&2
  repair_index
  sleep "0.$(( (RANDOM % 9) + 1 ))"
  attempt=$((attempt + 1))
done
