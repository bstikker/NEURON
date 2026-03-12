#!/usr/bin/env Rscript
# check_report_inputs.R
# Usage: Rscript check_report_inputs.R TVU11-7116

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript check_report_inputs.R <sample>")

sample <- args[1]

# Define all expected input files for the report
input_files <- c(
  paste0("results/", sample, "/merged_probes_methyl_calls_general.pdf"),
  paste0("results/", sample, "/MGMT_analysis/mgmt_prediction.tsv"),
  paste0("results/", sample, "/QDNAseq_ACE/", sample, "_segmented.png"),
  paste0("results/", sample, "/QDNAseq_ACE/ACE_summary.tsv"),
  paste0("results/", sample, "/pseudotime/confidence_vs_pseudotime.png"),
  paste0("results/", sample, "/pseudotime/confident_time.csv"),
  paste0("results/", sample, "/pseudotime/stable_time.csv")
)

# Check existence
missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  message("ERROR: The following input files are missing for sample ", sample, ":")
  for (f in missing_files) message("  ", f)
  stop("Cannot render report: missing input files.")
} else {
  message("All input files exist for sample ", sample)
}
