#!/usr/bin/env bash
# Install a standalone copy of the application and a PATH launcher.
set -euo pipefail
prefix="${HOME}/.local"
rscript="$(command -v Rscript || true)"
rc_files=()
usage() {
  cat <<'HELP'
Usage: bash scripts/install-cli.sh [options]
  --prefix DIR       Install under DIR/bin and DIR/lib/cytof-qc (default ~/.local)
  --rscript FILE     Rscript executable to use (default: command -v Rscript)
  --shell-rc FILE    Add an idempotent PATH entry to FILE; repeat for multiple files
  --help             Show this help
R dependencies are separate: run the selected Rscript with scripts/install.R.
HELP
}
while (($#)); do
  case "$1" in
    --help|-h) usage; exit 0 ;;
    --prefix|--rscript|--shell-rc)
      (($# >= 2)) && [[ -n "$2" ]] || { echo "Missing value for $1" >&2; exit 1; }
      case "$1" in
        --prefix) prefix="$2" ;;
        --rscript) rscript="$2" ;;
        --shell-rc) rc_files+=("$2") ;;
      esac
      shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done
[[ -n "$prefix" && "$prefix" != / ]] || { echo 'Choose a non-root prefix.' >&2; exit 1; }
[[ -x "$rscript" && ! -d "$rscript" ]] || { echo 'Rscript executable not found; use --rscript /path/to/Rscript.' >&2; exit 1; }
# Canonicalize directory components without resolving Rscript symlinks.
rscript="$(cd -- "$(dirname -- "$rscript")" && pwd -P)/$(basename -- "$rscript")"
"$rscript" --version
source_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
mkdir -p -- "$prefix"
prefix="$(cd -- "$prefix" && pwd -P)"
app="$prefix/lib/cytof-qc"
launcher="$prefix/bin/cytof-qc"
marker='# Managed by CyTOF QC CLI installer'
if [[ -e "$app" || -L "$app" ]]; then
  [[ ! -L "$app" && -f "$app/.cytof-qc-install" ]] || { echo "Refusing to replace unmanaged path: $app" >&2; exit 1; }
fi
if [[ -e "$launcher" || -L "$launcher" ]]; then
  [[ ! -L "$launcher" && -f "$launcher" ]] && grep -Fxq "$marker" "$launcher" || { echo "Refusing to replace unmanaged launcher: $launcher" >&2; exit 1; }
fi
mkdir -p -- "$prefix/bin" "$prefix/lib"
stage="$(mktemp -d "$prefix/lib/.cytof-qc.XXXXXX")"
trap 'rm -rf -- "$stage"' EXIT
mkdir "$stage/bin"
cp -R -- "$source_root/R" "$stage/R"
cp -- "$source_root/bin/cytof-qc.R" "$stage/bin/"
cp -- "$source_root/README.md" "$stage/README.md"
printf '%s\n' "$marker" > "$stage/.cytof-qc-install"
# Smoke test the copied application before replacing an installed version.
"$rscript" "$stage/bin/cytof-qc.R" --help >/dev/null
if [[ -d "$app" ]]; then rm -rf -- "$app"; fi
mv -- "$stage" "$app"
{
  printf '#!/usr/bin/env bash\n%s\n' "$marker"
  printf 'exec %q %q "$@"\n' "$rscript" "$app/bin/cytof-qc.R"
} > "$launcher"
chmod +x "$launcher"
# Single-quote paths for startup files shared by bash and POSIX login shells.
quoted_bin="'$(printf '%s' "$prefix/bin" | sed "s/'/'\\\\''/g")'"
path_line="export PATH=$quoted_bin:\"\$PATH\" # cytof-qc PATH"
for rc in "${rc_files[@]}"; do
  mkdir -p -- "$(dirname -- "$rc")"
  touch -- "$rc"
  if ! grep -Fxq -- "$path_line" "$rc"; then printf '\n%s\n' "$path_line" >> "$rc"; fi
done
printf '\nInstalled: %s\nR executable: %s\n' "$launcher" "$rscript"
printf 'For this shell, run:\n  export PATH=%s:"$PATH"\nThen run:\n  cytof-qc --help\n' "$quoted_bin"
