#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(png)
  library(scales)
})

# ------------------------ Args ------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript aggregate_sturgeon_timecourse_line_top10.R <sample_dir> <output_pdf> <summary_tsv> [read_counts_tsv]")
}

sample_dir      <- normalizePath(args[1], mustWork = TRUE)
output_pdf      <- args[2]
summary_tsv     <- args[3]
read_counts_tsv <- if (length(args) >= 4) args[4] else NA_character_
sample_id       <- basename(sample_dir)

safe_read_tsv <- function(f) {
  tryCatch(
    read.delim(f, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      warning(sprintf("Could not read %s: %s", f, e$message))
      NULL
    }
  )
}

infer_time_from_path <- function(path, sample_id) {
  m <- regexpr(paste0(sample_id, "_([0-9]+)min"), path, perl = TRUE)
  if (m[1] == -1) return(NA_real_)
  hit <- regmatches(path, m)[1]
  as.numeric(sub(paste0(sample_id, "_([0-9]+)min"), "\\1", hit, perl = TRUE))
}

# ------------------------ Load and reshape Sturgeon CSVs ------------------------
csv_files <- list.files(
  sample_dir,
  pattern = "merged_probes_methyl_calls_general.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(csv_files) == 0) {
  stop("No Sturgeon CSVs found in sample_dir")
}

parse_sturgeon_csv <- function(f, sample_id) {
  df <- tryCatch(
    read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      warning(sprintf("Could not read %s: %s", f, e$message))
      NULL
    }
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Expected current Sturgeon format:
  # one row, first col = number_probes, remaining cols = class scores
  if (!("number_probes" %in% colnames(df))) {
    warning(sprintf("Skipping %s because 'number_probes' column is missing", f))
    return(NULL)
  }

  time_min <- infer_time_from_path(f, sample_id)
  if (is.na(time_min)) {
    # skip the root-level non-timecourse CSV if present
    warning(sprintf("Skipping %s because timepoint could not be inferred from path", f))
    return(NULL)
  }

  cpg_used <- suppressWarnings(as.numeric(df$number_probes[1]))
  class_cols <- setdiff(colnames(df), "number_probes")
  if (length(class_cols) == 0) {
    warning(sprintf("Skipping %s because no class score columns were found", f))
    return(NULL)
  }

  scores <- suppressWarnings(as.numeric(df[1, class_cols]))
  names(scores) <- class_cols

  long_df <- data.frame(
    time_min = time_min,
    predicted_class = class_cols,
    confidence = scores,
    cpg_used = cpg_used,
    source_csv = basename(f),
    stringsAsFactors = FALSE
  )

  long_df
}

long_list <- lapply(csv_files, parse_sturgeon_csv, sample_id = sample_id)
long_list <- Filter(Negate(is.null), long_list)

if (length(long_list) == 0) {
  stop("No readable time-resolved Sturgeon CSVs found in sample_dir")
}

plot_df <- bind_rows(long_list) %>%
  mutate(
    time_min = as.numeric(time_min),
    confidence = as.numeric(confidence),
    predicted_class = as.character(predicted_class),
    cpg_used = as.numeric(cpg_used)
  ) %>%
  filter(!is.na(time_min), !is.na(confidence)) %>%
  arrange(time_min, desc(confidence))

if (nrow(plot_df) == 0) {
  stop("No valid Sturgeon rows remained after parsing")
}

# Per-timepoint summary: top class and top confidence
summary_table <- plot_df %>%
  group_by(time_min) %>%
  slice_max(order_by = confidence, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(time_min, confidence, predicted_class, cpg_used, source_csv)

# ------------------------ Merge read counts ------------------------
reads_df <- NULL
if (!is.na(read_counts_tsv) && file.exists(read_counts_tsv)) {
  reads_df <- safe_read_tsv(read_counts_tsv)
}

if (!is.null(reads_df) && all(c("time_min", "estimated_reads") %in% colnames(reads_df))) {
  reads_df <- reads_df %>%
    mutate(
      time_min = as.numeric(time_min),
      estimated_reads = as.numeric(estimated_reads)
    ) %>%
    select(time_min, estimated_reads) %>%
    distinct() %>%
    arrange(time_min)

  summary_table <- summary_table %>%
    left_join(reads_df, by = "time_min")
} else {
  summary_table$estimated_reads <- NA_real_
  reads_df <- data.frame(
    time_min = sort(unique(summary_table$time_min)),
    estimated_reads = NA_real_
  )
}

# ------------------------ Write summary TSV ------------------------
write.table(
  summary_table[, c("time_min", "estimated_reads", "confidence", "predicted_class", "cpg_used", "source_csv")],
  file = summary_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------ Locate CNV PNGs ------------------------
time_bins <- sort(unique(summary_table$time_min))

cnv_df <- data.frame(
  time_min = time_bins,
  cnv_png = sapply(time_bins, function(t) {
    file.path(
      sample_dir,
      paste0(sample_id, "_", t, "min"),
      "QDNAseq_ACE",
      paste0(sample_id, "_", t, "min_CNV.png")
    )
  }),
  stringsAsFactors = FALSE
) %>%
  mutate(exists = file.exists(cnv_png))

# ------------------------ Plot top 10 tumor classes ------------------------
top10_classes <- plot_df %>%
  group_by(predicted_class) %>%
  summarise(mean_conf = mean(confidence, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_conf)) %>%
  slice_head(n = 10) %>%
  pull(predicted_class)

plot_df_top10 <- plot_df %>%
  filter(predicted_class %in% top10_classes) %>%
  left_join(reads_df, by = "time_min") %>%
  arrange(time_min)

read_label_df <- plot_df_top10 %>%
  distinct(time_min, estimated_reads) %>%
  arrange(time_min)

p <- ggplot(plot_df_top10, aes(x = time_min, y = confidence, color = predicted_class)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  geom_hline(yintercept = 0.80, linetype = "dotted", color = "goldenrod2") +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "red") +
  geom_text(
    data = read_label_df,
    aes(x = time_min, y = 1.02, label = scales::comma(round(estimated_reads))),
    inherit.aes = FALSE,
    size = 3,
    angle = 45,
    hjust = 0,
    vjust = 0
  ) +
  annotate("text", x = min(read_label_df$time_min), y = 1.08,
           label = "Cumulative reads", hjust = 0, size = 3.5, fontface = "bold") +
  theme_minimal(base_size = 11) +
  labs(
    title = paste0(sample_id, ": Sturgeon confidence over time"),
    x = "Time (min)",
    y = "Sturgeon confidence",
    color = "Tumor class"
  ) +
  scale_y_continuous(limits = c(0, 1.1)) +
  scale_color_hue() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# ------------------------ Export PDF with CNV pages ------------------------
pdf(output_pdf, width = 10, height = 7)

print(p)

cnv_df_existing <- cnv_df %>%
  filter(exists) %>%
  arrange(time_min)

if (nrow(cnv_df_existing) > 0) {
  for (i in seq_len(nrow(cnv_df_existing))) {
    f <- cnv_df_existing$cnv_png[i]
    t <- cnv_df_existing$time_min[i]

    img <- tryCatch(readPNG(f), error = function(e) NULL)

    if (!is.null(img)) {
      grid.newpage()

      grid.text(
        paste0(sample_id, " - CNV profile at ", t, " min"),
        y = unit(0.97, "npc"),
        gp = gpar(fontsize = 14, fontface = "bold")
      )

      # Get page area to draw into (fractions of page)
      page_width_npc  <- 0.95  # fraction of page width
      page_height_npc <- 0.86  # fraction of page height
      page_aspect <- page_width_npc / page_height_npc

      img_h <- nrow(img)
      img_w <- ncol(img)
      img_aspect <- img_w / img_h

      # Scale while preserving aspect ratio
      if (img_aspect > page_aspect) {
        # image is wider than page box → width-limited
        draw_w <- page_width_npc
        draw_h <- page_width_npc / img_aspect
      } else {
        # image is taller than page box → height-limited
        draw_h <- page_height_npc
        draw_w <- page_height_npc * img_aspect
      }

      # Draw centered
      grid.raster(
        img,
        x = 0.5,
        y = 0.5,  # center vertically
        width = unit(draw_w, "npc"),
        height = unit(draw_h, "npc"),
        interpolate = FALSE
      )
    }
  }
}

dev.off()
cat("\nReport written to", output_pdf, "\n")
