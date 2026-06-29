#!/usr/bin/env bash
set -euo pipefail

FULL_BAM="${1:?Missing full BAM}"
READ_LIST="${2:?Missing read list}"
SAMPLE_ID="${3:?Missing sample ID}"
TIME_BIN="${4:?Missing time bin}"
PROJECT_ROOT="${5:-.}"

THREADS="${THREADS:-10}"

PROBES_FILE="${PROJECT_ROOT}/reference/probes/probelocs_chm13.bed"
MODEL_FILE="${PROJECT_ROOT}/reference/models/general.zip"

FINAL_DIR="${PROJECT_ROOT}/results/${SAMPLE_ID}/${SAMPLE_ID}_${TIME_BIN}min"
mkdir -p "${FINAL_DIR}"

BASE_TMP="${TMPDIR:-${PROJECT_ROOT}/tmp}"
mkdir -p "${BASE_TMP}"
WORKDIR="$(mktemp -d "${BASE_TMP}/${SAMPLE_ID}_${TIME_BIN}min_XXXXXX")"

cleanup() {
  rm -rf "${WORKDIR}"
}
trap cleanup EXIT

SUB_BAM="${WORKDIR}/${SAMPLE_ID}_${TIME_BIN}min.subset.bam"
ADJUSTED_BAM="${WORKDIR}/${SAMPLE_ID}_${TIME_BIN}min.adjusted.bam"
MODKIT_OUT="${WORKDIR}/modkit_extracted.txt"
STURGEON_IN="${WORKDIR}/sturgeon_input"
STURGEON_OUT="${WORKDIR}/sturgeon_output"

mkdir -p "${STURGEON_IN}" "${STURGEON_OUT}"

[[ -f "${FULL_BAM}" ]] || { echo "[ERROR] BAM not found: ${FULL_BAM}" >&2; exit 1; }
[[ -f "${READ_LIST}" ]] || { echo "[ERROR] Read list not found: ${READ_LIST}" >&2; exit 1; }
[[ -f "${PROBES_FILE}" ]] || { echo "[ERROR] Probes file not found: ${PROBES_FILE}" >&2; exit 1; }
[[ -f "${MODEL_FILE}" ]] || { echo "[ERROR] Model file not found: ${MODEL_FILE}" >&2; exit 1; }

echo "[INFO] Subsetting BAM..."
samtools view -@ "${THREADS}" -b -N "${READ_LIST}" -o "${SUB_BAM}" "${FULL_BAM}"
samtools index -@ "${THREADS}" "${SUB_BAM}"

echo "[INFO] Running modkit adjust-mods..."
modkit adjust-mods --convert h m "${SUB_BAM}" "${ADJUSTED_BAM}"

echo "[INFO] Running modkit extract..."
modkit extract full -t "${THREADS}" "${ADJUSTED_BAM}" "${MODKIT_OUT}"

cp "${MODKIT_OUT}" "${STURGEON_IN}/"

echo "[INFO] Running Sturgeon inputtobed..."
sturgeon inputtobed \
  --margin 50 \
  -i "${STURGEON_IN}" \
  -o "${STURGEON_OUT}" \
  -s modkit \
  --probes-file "${PROBES_FILE}"

echo "[INFO] Running Sturgeon predict..."
sturgeon predict \
  -p \
  --i "${STURGEON_OUT}" \
  -o "${STURGEON_OUT}" \
  -m "${MODEL_FILE}"

for f in merged_probes_methyl_calls_general.csv merged_probes_methyl_calls_general.pdf merged_probes_methyl_calls.bed; do
  if [[ -f "${STURGEON_OUT}/${f}" ]]; then
    cp -f "${STURGEON_OUT}/${f}" "${FINAL_DIR}/"
  fi
done

cat > "${FINAL_DIR}/run_metadata.tsv" <<EOF
sample_id	time_bin_min	full_bam	read_list
${SAMPLE_ID}	${TIME_BIN}	${FULL_BAM}	${READ_LIST}
EOF

echo "[INFO] Completed ${SAMPLE_ID} ${TIME_BIN} min"
