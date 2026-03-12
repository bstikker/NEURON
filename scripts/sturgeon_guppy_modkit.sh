#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <sturgeon_output_dir>" >&2
  exit 1
fi

OUTPUT_DIR="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PROBES_FILE="${REPO_ROOT}/reference/probes/probelocs_chm13.bed"
MODEL_FILE="${REPO_ROOT}/reference/models/general.zip"

if [[ ! -d "${OUTPUT_DIR}" ]]; then
  echo "Output directory not found: ${OUTPUT_DIR}" >&2
  exit 1
fi
if [[ ! -f "${PROBES_FILE}" ]]; then
  echo "Probes file not found: ${PROBES_FILE}" >&2
  exit 1
fi
if [[ ! -f "${MODEL_FILE}" ]]; then
  echo "Sturgeon model not found: ${MODEL_FILE}" >&2
  exit 1
fi

sturgeon inputtobed --margin 50 -i "${OUTPUT_DIR}" -o "${OUTPUT_DIR}" -s modkit --probes-file "${PROBES_FILE}"
sturgeon predict -p --i "${OUTPUT_DIR}" -o "${OUTPUT_DIR}" -m "${MODEL_FILE}"
