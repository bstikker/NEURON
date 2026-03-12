#!/usr/bin/env Rscript

# MGMT promoter methylation prediction (promoter-wide, with estimated coverage)
# Author: Bernard Stikker (modified by ChatGPT)
# Compatible with Snakemake
# --------------------------------------------------

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(readr))

# ---- Main Function ----
predict_mgmt_promoter <- function(modkit_file, model_file, promoter_bed,
                                  min_cpgs = 50, coverage_threshold = 3,
                                  verbose = TRUE) {
  
  # Load modkit methylation data
  modkit <- fread(modkit_file, select = c("chrom", "ref_position", "mod_qual"))

  # Load promoter BED
  promoter <- read_tsv(promoter_bed, col_names = c("chrom", "start", "end", "name"), show_col_types = FALSE) %>%
    select(chrom, start, end)

  # Filter modkit calls to promoter region
  modkit_filtered <- modkit %>%
    inner_join(promoter, by = "chrom") %>%
    filter(ref_position >= start & ref_position < end)

  # Handle edge case: no CpGs in region
  if (nrow(modkit_filtered) == 0) {
    return(list(
      avg_percent = NA,
      n_cpgs = 0,
      avg_coverage = NA,
      probability = NA,
      status = "Insufficient data",
      cpg_source = "promoter-wide",
      warning = "No CpG methylation calls found in promoter region"
    ))
  }

  # Compute CpG-level coverage estimate
  coverage_estimate <- modkit_filtered %>%
    group_by(chrom, ref_position) %>%
    summarise(cov = n(), .groups = "drop")

  avg_cov <- mean(coverage_estimate$cov, na.rm = TRUE)
  n_cpgs <- nrow(coverage_estimate)
  avg_beta <- mean(modkit_filtered$mod_qual, na.rm = TRUE)
  avg_percent <- avg_beta * 100

  # Warning messages
  warning_msgs <- c()
  if (n_cpgs < min_cpgs) {
    warning_msgs <- c(warning_msgs, sprintf("Low CpG count: %d (min %d)", n_cpgs, min_cpgs))
  }
  if (avg_cov < coverage_threshold) {
    warning_msgs <- c(warning_msgs, sprintf("Low average coverage: %.2fx", avg_cov))
  }
  warning_msgs <- c(warning_msgs, "Using promoter-wide CpGs; not limited to 137 model CpGs")

  # Predict using Patel model
  load(model_file)  # loads log.model
  probability <- predict(log.model, newdata = data.frame(average = avg_percent), type = "response")
  status <- ifelse(probability >= 0.5, "Methylated", "Unmethylated")

  # Verbose output
  if (verbose) {
    cat(sprintf("MGMT methylation: %.2f%% (n = %d CpGs, avg_cov = %.2fx)\n", avg_percent, n_cpgs, avg_cov))
    cat(sprintf("Prediction: %s (P = %.4f)\n", status, probability))
    if (length(warning_msgs) > 0) {
      cat("⚠️ Warnings:\n", paste0("- ", warning_msgs, collapse = "\n"), "\n")
    }
  }

  # Return result
  return(list(
    avg_percent = round(avg_percent, 2),
    n_cpgs = n_cpgs,
    avg_coverage = round(avg_cov, 2),
    probability = round(probability, 4),
    status = status,
    cpg_source = "promoter-wide",
    warning = paste(warning_msgs, collapse = "; ")
  ))
}

# ---- Argument Parser ----
option_list <- list(
  make_option(c("--modkit"), type = "character", help = "modkit-extracted methylation file"),
  make_option(c("--model"), type = "character", help = "Rdata file with logistic model"),
  make_option(c("--bed"), type = "character", help = "BED file for MGMT promoter"),
  make_option(c("--out"), type = "character", help = "Output result TSV"),
  make_option(c("--min_cpgs"), type = "integer", default = 50, help = "Minimum CpG count [default %default]"),
  make_option(c("--coverage_threshold"), type = "double", default = 3, help = "Minimum average coverage [default %default]"),
  make_option(c("--verbose"), action = "store_true", default = FALSE, help = "Verbose output")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ---- Run Analysis ----
result <- predict_mgmt_promoter(
  modkit_file = opt$modkit,
  model_file = opt$model,
  promoter_bed = opt$bed,
  min_cpgs = opt$min_cpgs,
  coverage_threshold = opt$coverage_threshold,
  verbose = opt$verbose
)

# ---- Write Result ----
write.table(as.data.frame(result), file = opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
