#!/usr/bin/env bash
set -euo pipefail
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
test_rscript="${CYTOF_TEST_RSCRIPT:-$(command -v Rscript)}"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
prefix="$tmp/prefix with space and ' quote"
rc="$tmp/bashrc"
bash "$root/scripts/install-cli.sh" --prefix "$prefix" --rscript "$test_rscript" --shell-rc "$rc" >/dev/null
(cd / && "$prefix/bin/cytof-qc" --help) | grep 'CyTOF QC' >/dev/null
bash "$root/scripts/install-cli.sh" --prefix "$prefix" --rscript "$test_rscript" --shell-rc "$rc" >/dev/null
[[ $(grep -c 'cytof-qc PATH' "$rc") == 1 ]]
bash -c 'source "$1"; cd /; cytof-qc --help' _ "$rc" | grep 'CyTOF QC' >/dev/null
[[ -f "$prefix/lib/cytof-qc/R/pipeline.R" ]]
if "$prefix/bin/cytof-qc" --unknown >/dev/null 2>&1; then echo 'Launcher lost failure status' >&2; exit 1; fi
mkdir -p "$tmp/conflict/bin"
printf 'unrelated\n' > "$tmp/conflict/bin/cytof-qc"
if bash "$root/scripts/install-cli.sh" --prefix "$tmp/conflict" --rscript "$test_rscript" >/dev/null 2>&1; then echo 'Overwrote unrelated command' >&2; exit 1; fi
[[ $(cat "$tmp/conflict/bin/cytof-qc") == unrelated ]]
echo 'PASS: CLI installation, independent working directory, quoting, reinstall, PATH and collision protection'
