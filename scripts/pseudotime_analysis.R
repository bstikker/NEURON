# pseudotime_analysis.R
# Usage: Rscript pseudotime_analysis.R <bam_path> <sample_id> <seq_summary>

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

if (length(args) < 3) {
  stop("Not enough arguments provided. Usage: Rscript pseudotime_analysis.R <bam_path> <sample_id> <seq_summary_file>")
}

bam_file <- normalizePath(args[1], mustWork = TRUE)
sample_id <- args[2]
seq_summary <- normalizePath(args[3], mustWork = TRUE)

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(data.table)
  library(stringr)
})

barcode_id <- stringr::str_extract(basename(bam_file), "barcode[0-9]+")
if (is.na(barcode_id)) stop("Could not extract barcode ID from BAM filename.")

cat("BAM:", bam_file, "\n")
cat("Sample ID:", sample_id, "\n")
cat("Summary file:", seq_summary, "\n")
cat("Extracted barcode ID:", barcode_id, "\n")
cat(sprintf("Running pseudotime analysis for %s\n", sample_id))

sturgeon_modkit_shellscript <- file.path(repo_root, "scripts", "sturgeon_guppy_modkit_SLURM.sh")
if (!file.exists(sturgeon_modkit_shellscript)) stop("Shell script not found: ", sturgeon_modkit_shellscript)

output_dir <- file.path(repo_root, "results", sample_id, "pseudotime")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seq_data <- read_tsv(seq_summary, show_col_types = FALSE) %>% arrange(start_time)
sturgeon_dir <- file.path(repo_root, "results", sample_id)
modkit_file <- file.path(sturgeon_dir, "modkit_extracted.txt")
if (!file.exists(modkit_file)) stop("modkit file not found: ", modkit_file)

modkit_data <- read_tsv(modkit_file, show_col_types = FALSE) %>%
  filter(ref_position != -1) %>%
  select(read_id, ref_position)

cpg_per_read <- modkit_data %>%
  group_by(read_id) %>%
  summarise(cpg_covered = n(), .groups = "drop") %>%
  arrange(match(read_id, seq_data$read_id))

non_cumulative_cpg <- c(51924, 104073, 124078, 149111, 173504, 194399, 207456, 217193, 232101, 241278, 247600, 258197)
pseudotime_bins <- data.frame(time_min = seq(5, 60, by = 5), expected_cpg = cumsum(non_cumulative_cpg))
cpg_counts <- data.frame(time_min = integer(), covered_cpg = integer())
confidence_all <- data.frame()
unique_cpg <- c()

for (i in seq_len(nrow(pseudotime_bins))) {
  bin_dir <- file.path(output_dir, paste0("bin_", pseudotime_bins$time_min[i]))
  dir.create(bin_dir, recursive = TRUE, showWarnings = FALSE)

  cumulative_reads <- character(0)
  cumulative_cpg <- 0
  for (row in seq_len(nrow(cpg_per_read))) {
    read_id <- cpg_per_read$read_id[row]
    cpg_count <- cpg_per_read$cpg_covered[row]
    cumulative_reads <- c(cumulative_reads, read_id)
    cumulative_cpg <- cumulative_cpg + cpg_count
    if (cumulative_cpg >= pseudotime_bins$expected_cpg[i]) break
  }

  reads_file <- file.path(bin_dir, paste0("reads_bin_", pseudotime_bins$time_min[i], ".txt"))
  writeLines(cumulative_reads, reads_file)

  bin_bam <- file.path(bin_dir, paste0("bin_", pseudotime_bins$time_min[i], ".bam"))
  adjusted_bam <- file.path(bin_dir, paste0("bin_", pseudotime_bins$time_min[i], "_adjusted.bam"))
  modkit_out <- file.path(bin_dir, paste0("bin_", pseudotime_bins$time_min[i], "_modkit.txt"))

  status <- system(sprintf(
    "samtools view -b -N %s %s -o %s",
    shQuote(reads_file), shQuote(bam_file), shQuote(bin_bam)
  ))
  if (status != 0) stop("samtools view failed for bin ", pseudotime_bins$time_min[i])

  status <- system(sprintf("samtools index %s", shQuote(bin_bam)))
  if (status != 0) stop("samtools index failed for bin ", pseudotime_bins$time_min[i])

  status <- system(sprintf(
    "modkit adjust-mods --convert h m %s %s",
    shQuote(bin_bam), shQuote(adjusted_bam)
  ))
  if (status != 0) stop("modkit adjust-mods failed for bin ", pseudotime_bins$time_min[i])

  status <- system(sprintf(
    "modkit extract full -t 20 %s %s",
    shQuote(adjusted_bam), shQuote(modkit_out)
  ))
  if (status != 0) stop("modkit extract failed for bin ", pseudotime_bins$time_min[i])

  if (file.exists(modkit_out)) {
    modkit_bin_data <- read_tsv(modkit_out, col_names = TRUE, show_col_types = FALSE)
    unique_cpg <- unique(c(unique_cpg, modkit_bin_data$ref_position))
  }

  cpg_counts <- rbind(
    cpg_counts,
    data.frame(time_min = pseudotime_bins$time_min[i], covered_cpg = length(unique_cpg))
  )

  status <- system(sprintf("%s %s", shQuote(sturgeon_modkit_shellscript), shQuote(bin_dir)))
  if (status != 0) stop("Sturgeon shell script failed for bin ", pseudotime_bins$time_min[i])

  confidence_file <- file.path(bin_dir, "merged_probes_methyl_calls_general.csv")
  if (file.exists(confidence_file)) {
    confidence_data <- read_csv(confidence_file, col_names = TRUE, show_col_types = FALSE)
    confidence_data_long <- pivot_longer(
      confidence_data,
      cols = -number_probes,
      names_to = "tumor_class",
      values_to = "confidence"
    )
    confidence_data_long <- confidence_data_long %>%
      mutate(
        time_min = pseudotime_bins$time_min[i],
        cpg_probes = confidence_data$number_probes[1],
        covered_cpg = cpg_counts$covered_cpg[cpg_counts$time_min == pseudotime_bins$time_min[i]]
      )
    confidence_all <- rbind(confidence_all, confidence_data_long)
  }
}

write_csv(confidence_all, file.path(output_dir, "confidence_vs_pseudotime.csv"))
final_point <- confidence_all %>%
  filter(time_min == max(time_min)) %>%
  slice_max(order_by = confidence, n = 1)

p <- ggplot(confidence_all, aes(x = cpg_probes, y = confidence, color = tumor_class)) +
  geom_line() +
  geom_point() +
  geom_text(data = final_point, aes(label = tumor_class), hjust = 1.2, vjust = 0, size = 5) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray") +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "gray") +
  labs(
    title = "Classification Confidence Over CpG probes covered",
    x = "CpG probes covered",
    y = "Confidence Score"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "confidence_vs_pseudotime.png"), p, width = 12, height = 10, bg = "white")
cat("Pseudotime analysis completed.\n")
