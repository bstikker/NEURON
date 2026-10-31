#!/usr/bin/env bash
# Collect software versions and reference/model checksums for the report.
# Runs inside the pipeline container. Writes a TSV to stdout with columns:
#   category  component  version  details
#
# Usage: collect_software_versions.sh <pipeline_version> <git_commit> <snakemake_version> [file ...]
# Every [file] is listed under "Models & references" with its MD5 checksum.
# A lookup that fails is reported as "not found"; it never stops the report.

set -uo pipefail

PIPELINE_VERSION="${1:-unknown}"
GIT_COMMIT="${2:-unknown}"
SNAKEMAKE_VERSION="${3:-unknown}"
shift 3 || true

row() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "${3:-not found}" "${4:-}"; }

# Version (and git commit, if installed from git) of a Python package in a given env.
py_pkg() {
  local python="$1" pkg="$2"
  [[ -x "$python" ]] || python="$(command -v python3 || command -v python)"
  "$python" - "$pkg" 2>/dev/null <<'EOF'
import json, sys
from importlib import metadata
pkg = sys.argv[1]
try:
    dist = metadata.distribution(pkg)
except metadata.PackageNotFoundError:
    sys.exit(1)
commit = ""
try:
    info = json.loads(dist.read_text("direct_url.json") or "{}")
    commit = info.get("vcs_info", {}).get("commit_id", "")
    if commit:
        commit = "commit " + commit[:7] + " (" + info.get("url", "") + ")"
except Exception:
    pass
print(dist.version + "\t" + commit)
EOF
}

printf 'category\tcomponent\tversion\tdetails\n'

# --- Pipeline ---
row "Pipeline" "NEURON" "$PIPELINE_VERSION" "git commit $GIT_COMMIT"
row "Pipeline" "Snakemake" "$SNAKEMAKE_VERSION" "workflow manager (runs outside the container)"

container_version=""
for labels in /.singularity.d/labels.json /.apptainer.d/labels.json; do
  if [[ -f "$labels" ]]; then
    container_version="$(sed -n 's/.*"Version"[^"]*"\([^"]*\)".*/\1/p' "$labels" | head -n 1)"
    break
  fi
done
container_file="${APPTAINER_CONTAINER:-${SINGULARITY_CONTAINER:-}}"
row "Pipeline" "Container" "$container_version" "$(basename "${container_file:-}" 2>/dev/null)"

# --- Classifiers and command-line tools ---
IFS=$'\t' read -r v d < <(py_pkg /opt/conda/envs/sturgeon_env/bin/python sturgeon || true)
row "Classifiers" "Sturgeon (v1)" "${v:-}" "${d:-}"
IFS=$'\t' read -r v d < <(py_pkg /opt/conda/envs/sturgeon_v2_env/bin/python sturgeon-v2 || true)
row "Classifiers" "Sturgeon v2" "${v:-}" "${d:-}"

v="$(modkit --version 2>/dev/null | awk '{print $NF}' | head -n 1)"
row "Tools" "modkit" "$v"
v="$(samtools --version 2>/dev/null | head -n 1 | awk '{print $2}')"
row "Tools" "samtools" "$v"
for env in sturgeon_env sturgeon_v2_env; do
  v="$(/opt/conda/envs/$env/bin/python --version 2>/dev/null | awk '{print $2}')"
  row "Tools" "Python ($env)" "$v"
done

# --- R and R packages ---
Rscript --vanilla -e '
pkgs <- c("QDNAseq", "CGHcall", "DNAcopy", "ACE", "Biobase", "GenomicRanges",
          "rmarkdown", "knitr", "kableExtra", "ggplot2", "dplyr", "readr", "data.table", "optparse")
cat(sprintf("R\tR\t%s.%s\t%s\n", R.version$major, R.version$minor, R.version$platform))
for (p in pkgs) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not found")
  cat(sprintf("R\t%s\t%s\t\n", p, v))
}' 2>/dev/null || row "R" "R" "" "Rscript failed"

v="$(xelatex --version 2>/dev/null | head -n 1)"
row "Report rendering" "XeTeX" "${v#XeTeX }"
v="$(pandoc --version 2>/dev/null | head -n 1 | awk '{print $2}')"
row "Report rendering" "Pandoc" "$v"

# --- Models and reference files ---
for f in "$@"; do
  if [[ -f "$f" ]]; then
    row "Models & references" "$(basename "$f")" "md5 $(md5sum "$f" | awk '{print $1}')" "$f"
  else
    row "Models & references" "$(basename "$f")" "not found" "$f"
  fi
done
