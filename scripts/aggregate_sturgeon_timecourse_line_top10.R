#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript aggregate_sturgeon_timecourse_line_top10.R <sample_dir> <output_pdf> <summary_tsv> [read_counts_tsv]")
}

sample_dir      <- normalizePath(args[1], mustWork = TRUE)
output_pdf      <- args[2]
summary_tsv     <- args[3]
read_counts_tsv <- if (length(args) >= 4) args[4] else NA_character_

cat("Sample directory:", sample_dir, "\n")

extract_minutes <- function(x) {
  x <- as.character(x)
  has_match <- grepl("_[0-9]+min", x)
  out <- sub(".*_([0-9]+)min.*", "\\1", x)
  out[!has_match] <- NA_character_
  as.numeric(out)
}

safe_read_csv <- function(f) {
  tryCatch(
    read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      warning(sprintf("Could not read %s: %s", f, e$message))
      NULL
    }
  )
}

safe_read_tsv <- function(f) {
  tryCatch(
    read.delim(f, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      warning(sprintf("Could not read %s: %s", f, e$message))
      NULL
    }
  )
}

shorten_class <- function(x, max_chars = 45) {
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

all_dirs <- list.dirs(sample_dir, recursive = FALSE, full.names = TRUE)

sturgeon_dirs <- all_dirs[
  grepl(paste0("^", basename(sample_dir), "_[0-9]+min$"), basename(all_dirs))
]

if (length(sturgeon_dirs) == 0) {
  stop("No Sturgeon time-bin folders found under sample directory.")
}

sturgeon_csvs <- unlist(lapply(sturgeon_dirs, function(d) {
  list.files(
    d,
    pattern = "^merged_probes_methyl_calls_general\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
}))

if (length(sturgeon_csvs) == 0) {
  stop("No merged_probes_methyl_calls_general.csv files found inside Sturgeon time-bin folders.")
}

summary_list <- list()
long_list <- list()

for (f in sturgeon_csvs) {
  df <- safe_read_csv(f)
  if (is.null(df) || nrow(df) == 0) next

  folder_name <- basename(dirname(f))
  if (!grepl("_[0-9]+min$", folder_name)) {
    path_parts <- strsplit(normalizePath(f), .Platform$file.sep)[[1]]
    hit <- grep("_[0-9]+min$", path_parts)
    if (length(hit) > 0) folder_name <- path_parts[max(hit)]
  }

  time_min <- extract_minutes(folder_name)
  if (is.na(time_min)) next
  if (!("number_probes" %in% colnames(df))) next

  class_cols <- setdiff(colnames(df), "number_probes")
  if (length(class_cols) == 0) next

  row_max <- apply(df[, class_cols, drop = FALSE], 1, max, na.rm = TRUE)
  best_row_idx <- which.max(row_max)
  best_row <- df[best_row_idx, , drop = FALSE]

  class_scores <- as.numeric(best_row[, class_cols, drop = TRUE])
  names(class_scores) <- class_cols

  best_class <- names(class_scores)[which.max(class_scores)]
  best_score <- max(class_scores, na.rm = TRUE)
  n_probes <- as.numeric(best_row$number_probes)

  summary_list[[length(summary_list) + 1]] <- data.frame(
    time_min = time_min,
    confidence = best_score,
    predicted_class = best_class,
    cpg_used = n_probes,
    source_csv = f,
    stringsAsFactors = FALSE
  )

  long_list[[length(long_list) + 1]] <- data.frame(
    time_min = time_min,
    tumor_class = names(class_scores),
    score = as.numeric(class_scores),
    cpg_used = n_probes,
    stringsAsFactors = FALSE
  )
}

if (length(summary_list) == 0) {
  stop("No usable Sturgeon CSVs could be parsed.")
}

summary_table <- bind_rows(summary_list) %>%
  distinct(time_min, .keep_all = TRUE) %>%
  arrange(time_min)

plot_df <- bind_rows(long_list) %>%
  arrange(time_min)

if (!is.na(read_counts_tsv) && file.exists(read_counts_tsv)) {
  reads_df <- safe_read_tsv(read_counts_tsv)
  if (!is.null(reads_df) && all(c("time_min", "estimated_reads") %in% colnames(reads_df))) {
    summary_table <- summary_table %>%
      left_join(reads_df %>% select(time_min, estimated_reads), by = "time_min")
  } else {
    summary_table$estimated_reads <- NA_integer_
  }
} else {
  summary_table$estimated_reads <- NA_integer_
}

write.table(
  summary_table[, c("time_min", "estimated_reads", "confidence", "predicted_class", "cpg_used", "source_csv")],
  file = summary_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

top_classes <- plot_df %>%
  group_by(tumor_class) %>%
  summarise(max_score = max(score, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_score)) %>%
  slice_head(n = 10) %>%
  pull(tumor_class)

plot_top <- plot_df %>%
  filter(tumor_class %in% top_classes) %>%
  mutate(tumor_class_short = shorten_class(tumor_class))

cpg_df <- summary_table %>%
  distinct(time_min, cpg_used, estimated_reads) %>%
  arrange(time_min)

p_timecourse <- ggplot(
  plot_top,
  aes(x = time_min, y = score, color = tumor_class_short, group = tumor_class_short)
) +
  geom_hline(yintercept = 0.80, linetype = "dotted", linewidth = 0.6, color = "#E6AC00") +
  geom_hline(yintercept = 0.95, linetype = "dotted", linewidth = 0.6, color = "#D73027") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  geom_text(
    data = cpg_df,
    aes(x = time_min, y = -0.04, label = cpg_used),
    inherit.aes = FALSE,
    angle = 90,
    vjust = 0.5,
    size = 2.5,
    color = "black"
  ) +
  annotate("text", x = min(cpg_df$time_min), y = 0.805, label = "0.80", hjust = -0.1, vjust = -0.4, size = 3, color = "#E6AC00") +
  annotate("text", x = min(cpg_df$time_min), y = 0.955, label = "0.95", hjust = -0.1, vjust = -0.4, size = 3, color = "#D73027") +
  scale_y_continuous(
    name = "Sturgeon class score",
    limits = c(-0.08, 1.02),
    expand = c(0.01, 0.01)
  ) +
  labs(
    x = "Time (minutes)",
    title = paste0(basename(sample_dir), " Sturgeon timecourse"),
    subtitle = "Top 10 class trajectories; numbers below x-axis = CpGs/probes used",
    color = "Tumor class"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

cnv_dirs <- all_dirs[grepl("^QDNAseq_ACE_[0-9]+min$", basename(all_dirs))]

cnv_minutes <- extract_minutes(basename(cnv_dirs))
cnv_dirs_hourly <- cnv_dirs[!is.na(cnv_minutes) & cnv_minutes %% 60 == 0]

if (length(cnv_dirs_hourly) == 0 && length(cnv_dirs) > 0) {
  cnv_dirs_hourly <- cnv_dirs
}

find_segmented_cnv_png <- function(d) {
  hits <- list.files(
    d,
    pattern = "_segmented\\.png$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  if (length(hits) > 0) return(hits[1])
  NA_character_
}

if (length(cnv_dirs_hourly) == 0) {
  cnv_png_files <- character(0)
  cnv_times <- numeric(0)
} else {
  cnv_png_files <- unname(sapply(cnv_dirs_hourly, find_segmented_cnv_png, USE.NAMES = TRUE))
  cnv_png_files <- as.character(cnv_png_files)
  cnv_png_files <- cnv_png_files[!is.na(cnv_png_files) & nzchar(cnv_png_files)]

  if (length(cnv_png_files) == 0) {
    cnv_times <- numeric(0)
  } else {
    cnv_times <- extract_minutes(basename(cnv_png_files))
    keep <- !is.na(cnv_times)
    cnv_png_files <- cnv_png_files[keep]
    cnv_times <- cnv_times[keep]

    ord <- order(cnv_times)
    cnv_png_files <- cnv_png_files[ord]
    cnv_times <- cnv_times[ord]
  }
}

pdf(output_pdf, width = 11, height = 8.5)
print(p_timecourse)

if (length(cnv_png_files) > 0) {
  if (requireNamespace("png", quietly = TRUE)) {
    for (i in seq_along(cnv_png_files)) {
      pf <- cnv_png_files[i]
      img <- tryCatch(
        png::readPNG(pf),
        error = function(e) NULL
      )

      if (!is.null(img)) {
        grid.newpage()
        grid.draw(rasterGrob(img, interpolate = FALSE))
        grid.text(
          paste0("CNV profile — ", cnv_times[i], " min"),
          x = unit(0.02, "npc"),
          y = unit(0.98, "npc"),
          just = c("left", "top"),
          gp = gpar(fontsize = 10, fontface = "bold")
        )
      }
    }
  }
}

dev.off()
cat("\nTimecourse report saved to ", output_pdf, "\n", sep = "")
