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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript aggregate_sturgeon_v2_timecourse.R <sample_dir> <output_pdf> <summary_tsv> [read_counts_tsv]")
}

sample_dir      <- normalizePath(args[1], mustWork = TRUE)
output_pdf      <- args[2]
summary_tsv     <- args[3]
read_counts_tsv <- if (length(args) >= 4) args[4] else NA_character_
sample_id       <- basename(sample_dir)

safe_read_tsv <- function(f) {
  tryCatch(read.delim(f, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
}

infer_time_from_path <- function(path, sample_id) {
  m <- regexpr(paste0(sample_id, "_([0-9]+)min"), path, perl = TRUE)
  if (m[1] == -1) return(NA_real_)
  hit <- regmatches(path, m)[1]
  as.numeric(sub(paste0(sample_id, "_([0-9]+)min"), "\\1", hit, perl = TRUE))
}

parse_sturgeon_csv <- function(f, sample_id) {
  df <- tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  if (!("number_probes" %in% colnames(df))) return(NULL)

  time_min <- infer_time_from_path(f, sample_id)
  if (is.na(time_min)) return(NULL)

  class_cols <- setdiff(colnames(df), "number_probes")
  if (length(class_cols) == 0) return(NULL)

  scores <- suppressWarnings(as.numeric(df[1, class_cols]))
  data.frame(
    classifier = "Sturgeon",
    time_min = time_min,
    predicted_class = class_cols,
    confidence = scores,
    cpg_used = suppressWarnings(as.numeric(df$number_probes[1])),
    source_csv = basename(f),
    stringsAsFactors = FALSE
  )
}

parse_sturgeon_v2_csv <- function(f, sample_id) {
  df <- tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)

  time_min <- infer_time_from_path(f, sample_id)
  if (is.na(time_min)) return(NULL)

  cn <- colnames(df)
  lower_cn <- tolower(cn)

  probe_idx <- which(lower_cn %in% c("probes", "number_probes", "cpg_used", "cpg_probes"))
  probe_val <- NA_real_
  if (length(probe_idx) > 0) {
    probe_val <- suppressWarnings(as.numeric(df[1, probe_idx[1]]))
  }

  # Remove obvious metadata columns and use the rest as class-confidence columns
  drop_cols <- unique(c(
    probe_idx,
    which(lower_cn %in% c("sample", "sample_id", "id", "barcode", "source_csv", "time_min"))
  ))
  class_idx <- setdiff(seq_along(cn), drop_cols)
  if (length(class_idx) == 0) return(NULL)

  scores <- suppressWarnings(as.numeric(df[1, class_idx]))
  keep <- !is.na(scores)
  if (!any(keep)) return(NULL)

  data.frame(
    classifier = "Sturgeon v2",
    time_min = time_min,
    predicted_class = cn[class_idx][keep],
    confidence = scores[keep],
    cpg_used = probe_val,
    source_csv = basename(f),
    stringsAsFactors = FALSE
  )
}

sturgeon_files <- list.files(sample_dir, pattern = "merged_probes_methyl_calls_general.csv$", recursive = TRUE, full.names = TRUE)
v2_files       <- list.files(sample_dir, pattern = "sturgeon_v2_outcome.csv$", recursive = TRUE, full.names = TRUE)

sturgeon_df <- bind_rows(Filter(Negate(is.null), lapply(sturgeon_files, parse_sturgeon_csv, sample_id = sample_id)))
v2_df       <- bind_rows(Filter(Negate(is.null), lapply(v2_files, parse_sturgeon_v2_csv, sample_id = sample_id)))

all_df <- bind_rows(sturgeon_df, v2_df) %>%
  mutate(
    time_min = as.numeric(time_min),
    confidence = as.numeric(confidence),
    predicted_class = as.character(predicted_class),
    cpg_used = as.numeric(cpg_used)
  ) %>%
  filter(!is.na(time_min), !is.na(confidence))

if (nrow(all_df) == 0) {
  stop("No valid Sturgeon/Sturgeon v2 timecourse rows found in sample_dir")
}

reads_df <- NULL
if (!is.na(read_counts_tsv) && file.exists(read_counts_tsv)) {
  reads_df <- safe_read_tsv(read_counts_tsv)
}

if (!is.null(reads_df) && all(c("time_min", "estimated_reads") %in% colnames(reads_df))) {
  reads_df <- reads_df %>%
    mutate(time_min = as.numeric(time_min), estimated_reads = as.numeric(estimated_reads)) %>%
    select(time_min, estimated_reads) %>%
    distinct() %>%
    arrange(time_min)
} else {
  reads_df <- data.frame(time_min = sort(unique(all_df$time_min)), estimated_reads = NA_real_)
}

summary_table <- all_df %>%
  group_by(classifier, time_min) %>%
  slice_max(order_by = confidence, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(reads_df, by = "time_min") %>%
  select(classifier, time_min, estimated_reads, confidence, predicted_class, cpg_used, source_csv)

write.table(
  summary_table,
  file = summary_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_classifier <- function(df, title_text, ylab_text) {
  if (nrow(df) == 0) return(NULL)

  top_classes <- df %>%
    group_by(predicted_class) %>%
    summarise(mean_conf = mean(confidence, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_conf)) %>%
    slice_head(n = 10) %>%
    pull(predicted_class)

  plot_df <- df %>%
    filter(predicted_class %in% top_classes) %>%
    left_join(reads_df, by = "time_min") %>%
    arrange(time_min)

  read_label_df <- plot_df %>%
    distinct(time_min, estimated_reads) %>%
    arrange(time_min)

  ggplot(plot_df, aes(x = time_min, y = confidence, color = predicted_class)) +
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
      title = paste0(sample_id, ": ", title_text),
      x = "Time (min)",
      y = ylab_text,
      color = "Tumor class"
    ) +
    scale_y_continuous(limits = c(0, 1.1)) +
    scale_color_hue() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

cnv_files <- list.files(sample_dir, pattern = "_CNV\\.png$", recursive = TRUE, full.names = TRUE)
cnv_df <- data.frame(
  cnv_png = cnv_files,
  time_min = sapply(cnv_files, infer_time_from_path, sample_id = sample_id),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(time_min)) %>%
  arrange(time_min)

pdf(output_pdf, width = 10, height = 7)

p1 <- plot_classifier(filter(all_df, classifier == "Sturgeon"), "Sturgeon confidence over time", "Sturgeon confidence")
if (!is.null(p1)) print(p1)

p2 <- plot_classifier(filter(all_df, classifier == "Sturgeon v2"), "Sturgeon v2 confidence over time", "Sturgeon v2 confidence")
if (!is.null(p2)) print(p2)

if (nrow(cnv_df) > 0) {
  for (i in seq_len(nrow(cnv_df))) {
    f <- cnv_df$cnv_png[i]
    t <- cnv_df$time_min[i]
    img <- tryCatch(readPNG(f), error = function(e) NULL)

    if (!is.null(img)) {
      grid.newpage()
      grid.text(
        paste0(sample_id, " - CNV profile at ", t, " min"),
        y = unit(0.97, "npc"),
        gp = gpar(fontsize = 14, fontface = "bold")
      )

      page_width_npc  <- 0.95
      page_height_npc <- 0.86
      page_aspect <- page_width_npc / page_height_npc

      img_h <- nrow(img)
      img_w <- ncol(img)
      img_aspect <- img_w / img_h

      if (img_aspect > page_aspect) {
        draw_w <- page_width_npc
        draw_h <- page_width_npc / img_aspect
      } else {
        draw_h <- page_height_npc
        draw_w <- page_height_npc * img_aspect
      }

      grid.raster(
        img,
        x = 0.5,
        y = 0.5,
        width = unit(draw_w, "npc"),
        height = unit(draw_h, "npc"),
        interpolate = FALSE
      )
    }
  }
}

dev.off()
cat("\nReport written to", output_pdf, "\n")
