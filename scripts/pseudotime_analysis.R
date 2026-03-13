#!/usr/bin/env Rscript

# pseudotime_analysis.R
# Usage:
#   Rscript pseudotime_analysis.R <bam_path> <sample_id> <seq_summary> [barcode]
#
# Optional environment variables:
#   PSEUDOTIME_BARCODE            Explicit barcode override, e.g. "barcode01"
#   STURGEON_HELPER              Path to sturgeon helper shell script
#   PSEUDOTIME_TIMEPOINTS        Comma-separated minutes, e.g. "5,10,15,20,30,45,60"
#   PSEUDOTIME_NONCUM_CPG        Comma-separated non-cumulative CpG targets
#   PSEUDOTIME_THREADS           Threads for modkit extract (default: 20)
#
# Outputs:
#   results/<sample_id>/pseudotime/
#     - confidence_vs_pseudotime.csv
#     - confidence_vs_pseudotime.png
#     - pseudotime_metadata.tsv

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript pseudotime_analysis.R <bam_path> <sample_id> <seq_summary> [barcode]")
}

bam_path    <- normalizePath(args[1], mustWork = TRUE)
sample_id   <- args[2]
seqsum_path <- normalizePath(args[3], mustWork = TRUE)

barcode_arg <- if (length(args) >= 4 && nzchar(args[4])) args[4] else Sys.getenv("PSEUDOTIME_BARCODE", "")
barcode_arg <- trimws(barcode_arg)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(readr)
})

# ----------------------------- helpers --------------------------------------

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")

stopf <- function(...) stop(sprintf(...), call. = FALSE)

require_file <- function(path, label = NULL) {
  if (!file.exists(path)) {
    stopf("%s not found: %s", ifelse(is.null(label), "File", label), path)
  }
  normalizePath(path, mustWork = TRUE)
}

find_first_existing <- function(paths, label) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) stopf("%s not found. Checked:\n- %s", label, paste(paths, collapse = "\n- "))
  normalizePath(hit, mustWork = TRUE)
}

run_cmd <- function(cmd, args = character(), fail_msg = NULL) {
  msg("Running: ", paste(c(cmd, args), collapse = " "))
  out <- system2(cmd, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  if (status != 0L) {
    if (length(out)) message(paste(out, collapse = "\n"))
    if (is.null(fail_msg)) fail_msg <- paste("Command failed:", cmd)
    stopf("%s (exit status %s)", fail_msg, status)
  }
  invisible(out)
}

get_barcode <- function(bam_path, barcode_arg = "") {
  if (nzchar(barcode_arg)) {
    msg("Using explicit barcode: ", barcode_arg)
    return(barcode_arg)
  }

  bam_base <- basename(bam_path)
  barcode_id <- stringr::str_extract(bam_base, "barcode\\d+")

  if (is.na(barcode_id) || is.null(barcode_id)) {
    msg("No barcode detected in BAM filename. Proceeding without barcode filtering.")
    return(NA_character_)
  }

  msg("Detected barcode from BAM filename: ", barcode_id)
  barcode_id
}

detect_time_column <- function(dt) {
  candidates <- c("start_time", "start_time_s", "start_time_sec", "start_time_seconds")
  hit <- candidates[candidates %in% names(dt)][1]
  if (is.na(hit)) {
    stopf("Could not find a sequencing-summary time column. Expected one of: %s",
          paste(candidates, collapse = ", "))
  }
  hit
}

detect_probe_column <- function(dt) {
  preferred <- c("number_probes", "cpg_probes")
  hit <- preferred[preferred %in% names(dt)][1]
  if (!is.na(hit)) return(hit)

  numeric_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
  if (length(numeric_cols) == 0) {
    stopf("Could not determine the probe-count column in Sturgeon output.")
  }
  numeric_cols[1]
}

parse_int_vector_env <- function(var, default) {
  x <- Sys.getenv(var, "")
  if (!nzchar(x)) return(default)
  vals <- suppressWarnings(as.integer(trimws(strsplit(x, ",")[[1]])))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(default)
  vals
}

# ----------------------------- paths ----------------------------------------

bam_path <- require_file(bam_path, "BAM")
seqsum_path <- require_file(seqsum_path, "Sequencing summary")

sturgeon_helper <- Sys.getenv("STURGEON_HELPER", "")
if (!nzchar(sturgeon_helper)) {
  sturgeon_helper <- find_first_existing(
    c(
      file.path(repo_root, "scripts", "sturgeon_guppy_modkit_SLURM.sh"),
      file.path(repo_root, "scripts", "sturgeon_guppy_modkit.sh")
    ),
    "Sturgeon helper script"
  )
} else {
  sturgeon_helper <- require_file(sturgeon_helper, "Sturgeon helper script")
}

sturgeon_dir <- file.path(repo_root, "results", sample_id)
output_dir <- file.path(sturgeon_dir, "pseudotime")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

modkit_file <- file.path(sturgeon_dir, "modkit_extracted.txt")
modkit_file <- require_file(modkit_file, "modkit extracted file")

# ----------------------------- config ---------------------------------------

barcode_id <- get_barcode(bam_path, barcode_arg)

default_timepoints <- seq(5L, 60L, by = 5L)
default_noncum_cpg <- c(51924L, 104073L, 124078L, 149111L, 173504L, 194399L,
                        207456L, 217193L, 232101L, 241278L, 247600L, 258197L)

timepoints_min <- parse_int_vector_env("PSEUDOTIME_TIMEPOINTS", default_timepoints)
noncum_cpg     <- parse_int_vector_env("PSEUDOTIME_NONCUM_CPG", default_noncum_cpg)

if (length(timepoints_min) != length(noncum_cpg)) {
  stopf("Length mismatch: PSEUDOTIME_TIMEPOINTS has %d values, PSEUDOTIME_NONCUM_CPG has %d values.",
        length(timepoints_min), length(noncum_cpg))
}

modkit_threads <- suppressWarnings(as.integer(Sys.getenv("PSEUDOTIME_THREADS", "20")))
if (is.na(modkit_threads) || modkit_threads < 1L) modkit_threads <- 20L

pseudotime_bins <- data.frame(
  time_min = timepoints_min,
  expected_cpg = cumsum(noncum_cpg)
)

msg("BAM: ", bam_path)
msg("Sample ID: ", sample_id)
msg("Sequencing summary: ", seqsum_path)
msg("Results dir: ", output_dir)
msg("Using helper script: ", sturgeon_helper)

# ----------------------- read sequencing summary ----------------------------

seq_data <- data.table::fread(seqsum_path, showProgress = FALSE)
if (!"read_id" %in% names(seq_data)) {
  stopf("Sequencing summary lacks required column: read_id")
}

time_col <- detect_time_column(seq_data)

if (!is.na(barcode_id) && "barcode_arrangement" %in% names(seq_data)) {
  n_before <- nrow(seq_data)
  seq_data <- seq_data[barcode_arrangement == barcode_id]
  msg("Barcode filtering kept ", nrow(seq_data), " / ", n_before, " rows from sequencing summary.")
} else {
  msg("Skipping barcode filtering in sequencing summary.")
}

if (nrow(seq_data) == 0) {
  stopf("No sequencing-summary rows remain after filtering.")
}

setorderv(seq_data, time_col)

# ----------------------- read modkit extracted file -------------------------

modkit_data <- data.table::fread(modkit_file, showProgress = FALSE)

required_modkit_cols <- c("read_id", "ref_position")
missing_modkit_cols <- setdiff(required_modkit_cols, names(modkit_data))
if (length(missing_modkit_cols) > 0) {
  stopf("modkit extracted file is missing required columns: %s",
        paste(missing_modkit_cols, collapse = ", "))
}

modkit_data <- modkit_data[ref_position != -1]

if (nrow(modkit_data) == 0) {
  stopf("No valid rows in modkit extracted file after filtering ref_position != -1.")
}

# Count CpG observations per read from the full-sample modkit output
cpg_per_read <- modkit_data[, .(cpg_covered = .N), by = read_id]

# Keep only reads that appear in sequencing summary
cpg_per_read <- cpg_per_read[read_id %in% seq_data$read_id]

if (nrow(cpg_per_read) == 0) {
  stopf("No overlap between modkit read IDs and sequencing-summary read IDs.")
}

# Preserve sequencing-summary order
seq_order <- data.table(read_id = seq_data$read_id, seq_idx = seq_len(nrow(seq_data)))
cpg_per_read <- merge(cpg_per_read, seq_order, by = "read_id", all.x = TRUE)
cpg_per_read <- cpg_per_read[!is.na(seq_idx)]
setorder(cpg_per_read, seq_idx)

# Cumulative CpG target selector
cpg_per_read[, cum_cpg := cumsum(cpg_covered)]

# ----------------------------- iterate bins ---------------------------------

confidence_all <- list()
metadata_rows  <- list()
covered_cpg_tracker <- integer(0)

for (i in seq_len(nrow(pseudotime_bins))) {
  tp_min <- pseudotime_bins$time_min[i]
  target_cpg <- pseudotime_bins$expected_cpg[i]

  msg("Processing pseudotime bin ", tp_min, " min (target cumulative CpG = ", target_cpg, ")")

  bin_dir <- file.path(output_dir, paste0("bin_", tp_min))
  dir.create(bin_dir, recursive = TRUE, showWarnings = FALSE)

  # Select reads up to this cumulative CpG target
  cutoff_idx <- which(cpg_per_read$cum_cpg >= target_cpg)[1]
  if (is.na(cutoff_idx)) {
    cutoff_idx <- nrow(cpg_per_read)
    msg("Target exceeds available cumulative CpG. Using all remaining reads.")
  }

  selected_reads <- cpg_per_read$read_id[seq_len(cutoff_idx)]
  selected_cum_cpg <- cpg_per_read$cum_cpg[cutoff_idx]

  if (length(selected_reads) == 0) {
    msg("No reads selected for bin ", tp_min, "; skipping.")
    next
  }

  reads_file   <- file.path(bin_dir, sprintf("reads_bin_%s.txt", tp_min))
  bin_bam      <- file.path(bin_dir, sprintf("bin_%s.bam", tp_min))
  adjusted_bam <- file.path(bin_dir, sprintf("bin_%s_adjusted.bam", tp_min))
  modkit_out   <- file.path(bin_dir, sprintf("bin_%s_modkit.txt", tp_min))
  conf_file    <- file.path(bin_dir, "merged_probes_methyl_calls_general.csv")

  writeLines(selected_reads, reads_file)

  run_cmd(
    "samtools",
    c("view", "-b", "-N", reads_file, bam_path, "-o", bin_bam),
    fail_msg = sprintf("samtools view failed for pseudotime bin %s", tp_min)
  )

  run_cmd(
    "samtools",
    c("index", bin_bam),
    fail_msg = sprintf("samtools index failed for pseudotime bin %s", tp_min)
  )

  run_cmd(
    "modkit",
    c("adjust-mods", "--convert", "h", "m", bin_bam, adjusted_bam),
    fail_msg = sprintf("modkit adjust-mods failed for pseudotime bin %s", tp_min)
  )

  run_cmd(
    "modkit",
    c("extract", "full", "-t", as.character(modkit_threads), adjusted_bam, modkit_out),
    fail_msg = sprintf("modkit extract failed for pseudotime bin %s", tp_min)
  )

  if (file.exists(modkit_out)) {
    modkit_bin_data <- data.table::fread(modkit_out, showProgress = FALSE)
    if ("ref_position" %in% names(modkit_bin_data)) {
      covered_cpg_tracker <- unique(c(covered_cpg_tracker, modkit_bin_data$ref_position))
    }
  }

  run_cmd(
    sturgeon_helper,
    c(bin_dir),
    fail_msg = sprintf("Sturgeon helper script failed for pseudotime bin %s", tp_min)
  )

  if (!file.exists(conf_file)) {
    msg("Sturgeon output not found for bin ", tp_min, ": ", conf_file)
    next
  }

  conf_dt <- data.table::fread(conf_file, showProgress = FALSE)
  if (nrow(conf_dt) == 0) {
    msg("Empty Sturgeon output for bin ", tp_min, "; skipping.")
    next
  }

  probe_col <- detect_probe_column(conf_dt)

  conf_long <- as.data.frame(conf_dt) |>
    tidyr::pivot_longer(
      cols = -all_of(probe_col),
      names_to = "tumor_class",
      values_to = "confidence"
    ) |>
    dplyr::mutate(
      time_min = tp_min,
      sample_id = sample_id,
      barcode = ifelse(is.na(barcode_id), "", barcode_id),
      cpg_probes = .data[[probe_col]][1],
      covered_cpg = length(covered_cpg_tracker),
      selected_reads = length(selected_reads),
      selected_cumulative_cpg = selected_cum_cpg
    ) |>
    dplyr::select(sample_id, barcode, time_min, cpg_probes, covered_cpg,
                  selected_reads, selected_cumulative_cpg, tumor_class, confidence)

  confidence_all[[length(confidence_all) + 1L]] <- conf_long

  best_idx <- which.max(conf_long$confidence)
  metadata_rows[[length(metadata_rows) + 1L]] <- data.frame(
    sample_id = sample_id,
    barcode = ifelse(is.na(barcode_id), "", barcode_id),
    time_min = tp_min,
    selected_reads = length(selected_reads),
    selected_cumulative_cpg = selected_cum_cpg,
    covered_cpg = length(covered_cpg_tracker),
    top_class = conf_long$tumor_class[best_idx],
    top_confidence = conf_long$confidence[best_idx],
    stringsAsFactors = FALSE
  )
}

if (length(confidence_all) == 0) {
  stopf("No pseudotime confidence data were generated.")
}

confidence_all <- bind_rows(confidence_all)
metadata_df <- bind_rows(metadata_rows)

out_csv <- file.path(output_dir, "confidence_vs_pseudotime.csv")
out_png <- file.path(output_dir, "confidence_vs_pseudotime.png")
out_meta <- file.path(output_dir, "pseudotime_metadata.tsv")

readr::write_csv(confidence_all, out_csv)
data.table::fwrite(metadata_df, out_meta, sep = "\t")

final_point <- confidence_all |>
  dplyr::filter(time_min == max(time_min, na.rm = TRUE)) |>
  dplyr::slice_max(order_by = confidence, n = 1, with_ties = FALSE)

p <- ggplot(confidence_all, aes(x = cpg_probes, y = confidence, color = tumor_class, group = tumor_class)) +
  geom_line() +
  geom_point() +
  geom_text(data = final_point, aes(label = tumor_class), hjust = 1.1, vjust = 0, size = 4.5, show.legend = FALSE) +
  geom_hline(yintercept = 0.80, linetype = "dotted", color = "gray50") +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "gray50") +
  labs(
    title = sprintf("Sturgeon classification confidence over pseudotime: %s", sample_id),
    x = "CpG probes covered (Sturgeon output)",
    y = "Confidence score"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(out_png, p, width = 12, height = 10, bg = "white")

msg("Pseudotime analysis completed.")
msg("CSV: ", out_csv)
msg("PNG: ", out_png)
msg("Metadata: ", out_meta)