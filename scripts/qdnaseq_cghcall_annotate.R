#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript scripts/qdnaseq_cghcall_annotate.R <sample_id> <segmented_rds> <gene_bed> <output_dir> [generate_locuszoom]")
}

sample_id          <- args[1]
segmented_rds      <- normalizePath(args[2], mustWork = TRUE)
gene_bed_file      <- normalizePath(args[3], mustWork = TRUE)
output_dir         <- args[4]
generate_locuszoom <- if (length(args) >= 5) tolower(args[5]) %in% c("true", "1", "yes") else FALSE

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(GenomicRanges)
  library(CGHcall)
})

first  <- dplyr::first
slice  <- dplyr::slice
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
summarise <- dplyr::summarise
rename <- dplyr::rename
transmute <- dplyr::transmute

safe_ggsave <- function(filename, plot, width = 6, height = 4, dpi = 200) {
  tmp <- tempfile(fileext = ".png")
  ggsave(tmp, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(tmp, filename, overwrite = TRUE)
  if (!ok) warning("Failed to write output file: ", filename)
  invisible(filename)
}

plot_gene_region_qc <- function(cn, genes_df, gene_symbol, sample_id, flank_mb = 10) {
  g <- genes_df %>% filter(name == gene_symbol) %>% slice(1)
  if (nrow(g) == 0) return(ggplot() + ggtitle(paste("Gene not found:", gene_symbol)))

  gene_chr <- g$chrom
  g_start  <- g$start
  g_end    <- g$end

  flank_bp   <- flank_mb * 1e6
  region_min <- max(1, g_start - flank_bp)
  region_max <- g_end + flank_bp

  fd <- fData(cn)
  bins_chr <- gsub("^chr", "", as.character(fd$chromosome))

  copynum_mat <- assayDataElement(cn, "copynumber")
  raw_copynum <- copynum_mat[, sampleNames(cn)[1]]
  log2ratio   <- log2(raw_copynum)
  log2ratio[!is.finite(log2ratio)] <- NA

  region_bins <- fd %>%
    mutate(chrom = bins_chr, log2ratio = log2ratio) %>%
    filter(chrom == gene_chr, start <= region_max, end >= region_min, !is.na(log2ratio)) %>%
    mutate(mid_pos = (start + end) / 2,
           overlaps_gene = (end >= g_start & start <= g_end))

  if (nrow(region_bins) == 0) return(ggplot() + ggtitle(paste("No bins for", gene_symbol)))

  median_shift <- median(region_bins$log2ratio, na.rm = TRUE)
  region_bins$log2ratio <- region_bins$log2ratio - median_shift

  ggplot(region_bins, aes(x = mid_pos / 1e6, y = log2ratio)) +
    geom_point(alpha = 0.5) +
    geom_point(data = subset(region_bins, overlaps_gene), color = "red", size = 2.5) +
    annotate("rect", xmin = g_start / 1e6, xmax = g_end / 1e6,
             ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1) +
    annotate("text", x = (g_start + g_end) / 2 / 1e6, y = Inf,
             label = gene_symbol, vjust = 1.5, fontface = "bold", size = 3) +
    xlab(paste0("Chr ", gene_chr, " (Mb)")) +
    ylab("log2 ratio") +
    ggtitle(paste0("Region: ", gene_symbol, " (", sample_id, ")")) +
    theme_bw()
}

message("Reading segmented RDS: ", segmented_rds)
cn_raw <- readRDS(segmented_rds)
message("Running CGHcall")
cn <- callBins(cn_raw, method = "CGHcall", organism = "human")
sample_name_in_obj <- sampleNames(cn)[1]

gene_ann_raw <- read_tsv(gene_bed_file, col_names = FALSE, show_col_types = FALSE)
colnames(gene_ann_raw) <- c("chrom_raw", "start", "end", "name", "binnr", "yoffset")

genes_df <- gene_ann_raw %>%
  mutate(
    chrom   = gsub("^chr", "", chrom_raw),
    start   = as.integer(start),
    end     = as.integer(end),
    binnr   = as.integer(binnr),
    yoffset = suppressWarnings(as.numeric(yoffset))
  ) %>%
  filter(!is.na(start), !is.na(end), !is.na(binnr))

calls_vec <- assayDataElement(cn, "calls")[, sample_name_in_obj]
probloss  <- assayDataElement(cn, "probloss")[, sample_name_in_obj]
probnorm  <- assayDataElement(cn, "probnorm")[, sample_name_in_obj]
probgain  <- assayDataElement(cn, "probgain")[, sample_name_in_obj]
probamp   <- assayDataElement(cn, "probamp")[, sample_name_in_obj]

call_prob <- ifelse(is.na(calls_vec), NA_real_,
                    ifelse(calls_vec <= -1, probloss,
                           ifelse(calls_vec == 0, probnorm,
                                  ifelse(calls_vec == 1, probgain,
                                         ifelse(calls_vec >= 2, probamp, NA_real_)))))
call_prob[is.na(call_prob)] <- 0.5

bins_all <- fData(cn) %>% as.data.frame()
copynum_mat <- assayDataElement(cn, "copynumber")
seg_mat     <- assayDataElement(cn, "segmented")

bins_all$log2ratio <- log2(copynum_mat[, sample_name_in_obj])
bins_all$segmented <- log2(seg_mat[, sample_name_in_obj])
bins_all$log2ratio[!is.finite(bins_all$log2ratio)] <- NA
bins_all$segmented[!is.finite(bins_all$segmented)] <- NA
bins_all$bin_global <- seq_len(nrow(bins_all))
bins_all$chromosome <- suppressWarnings(as.integer(as.character(bins_all$chromosome)))
bins_all$call <- as.integer(calls_vec)
bins_all$call_prob <- as.numeric(call_prob)

autosomes <- bins_all %>%
  filter(chromosome %in% 1:22, !is.na(log2ratio)) %>%
  arrange(chromosome, start) %>%
  mutate(bin_index = row_number()) %>%
  mutate(call_class = case_when(
    call <= -1 ~ "loss",
    call == 1  ~ "gain",
    call >= 2  ~ "amp",
    TRUE       ~ "neutral"
  ))

cn_used <- bins_all %>%
  filter(chromosome %in% 1:22, use == TRUE, !is.na(log2ratio)) %>%
  arrange(chromosome, start) %>%
  mutate(bin_index = match(bin_global, autosomes$bin_global))

gene_ann <- genes_df %>%
  mutate(bin_global = binnr,
         chromosome = as.integer(chrom))

autosomes <- autosomes %>%
  left_join(gene_ann %>% select(bin_global, name, yoffset), by = "bin_global") %>%
  mutate(is_gene_bin = !is.na(name))

chr_bounds <- autosomes %>%
  group_by(chromosome) %>%
  summarise(min_bin = min(bin_index), max_bin = max(bin_index), .groups = "drop")

chr_labels <- chr_bounds %>% mutate(label_pos = (min_bin + max_bin) / 2)

seg_bins <- cn_used %>%
  filter(!is.na(segmented)) %>%
  arrange(bin_global) %>%
  mutate(seg_group = cumsum(c(1, diff(segmented) != 0)))

segments_df <- seg_bins %>%
  group_by(seg_group) %>%
  summarise(
    chrom        = first(chromosome),
    start        = min(start),
    end          = max(end),
    n_bins       = n(),
    seg_value_raw= first(segmented),
    x_start      = min(bin_index),
    x_end        = max(bin_index),
    seg_call     = median(call, na.rm = TRUE),
    seg_prob     = median(call_prob, na.rm = TRUE),
    .groups      = "drop"
  ) %>%
  mutate(call_class = case_when(
    seg_call <= -1 ~ "loss",
    seg_call == 1  ~ "gain",
    seg_call >= 2  ~ "amp",
    TRUE           ~ "neutral"
  ))

pd <- pData(cn)
obs_sigma <- if ("expected.variance" %in% colnames(pd)) sqrt(pd[sample_name_in_obj, "expected.variance"]) else NA_real_
condition <- QDNAseq:::binsToUse(cn)
copynumber_cond <- assayDataElement(cn, "copynumber")[condition, sample_name_in_obj]
delta_sigma <- QDNAseq:::sdDiffTrim(copynumber_cond, na.rm = TRUE)
n_segments <- dplyr::n_distinct(seg_bins$seg_group)

n_bins_used   <- nrow(cn_used)
bin_width_bp  <- median(cn_used$end - cn_used$start, na.rm = TRUE)
bin_width_kbp <- round(bin_width_bp / 1000)

left_label <- sprintf("%dk x %d kbp, %d segments",
                      round(n_bins_used / 1000),
                      bin_width_kbp,
                      n_segments)
right_label <- sprintf("E σ = %.3f, Δσ̂ = %.3f", obs_sigma, delta_sigma)

available_assays <- assayDataElementNames(cn)
if ("total.reads" %in% colnames(pd)) {
  total_reads <- pd[sample_name_in_obj, "total.reads"]
} else if ("reads" %in% colnames(pd)) {
  total_reads <- pd[sample_name_in_obj, "reads"]
} else if ("counts" %in% available_assays) {
  total_reads <- sum(assayDataElement(cn, "counts")[, sample_name_in_obj], na.rm = TRUE)
} else {
  total_reads <- NA_real_
}
formatted_reads <- if (is.na(total_reads)) "reads" else paste0(format(total_reads, big.mark = ",", scientific = FALSE), " reads")
plot_title <- paste0(sample_id, " (", formatted_reads, ")")

median_shift <- median(cn_used$log2ratio, na.rm = TRUE)
autosomes$log2ratio <- autosomes$log2ratio - median_shift
segments_df$seg_value <- segments_df$seg_value_raw - median_shift

gene_bins <- autosomes %>%
  filter(is_gene_bin) %>%
  distinct(bin_index, name, log2ratio, yoffset)

x_min <- min(autosomes$bin_index, na.rm = TRUE)
x_max <- max(autosomes$bin_index, na.rm = TRUE)

p <- ggplot(autosomes, aes(x = bin_index, y = log2ratio)) +
  geom_point(aes(color = call_class, alpha = call_prob), size = 0.6) +
  geom_point(data = gene_bins, aes(x = bin_index, y = log2ratio),
             inherit.aes = FALSE, color = "red", size = 1.2) +
  geom_text(data = gene_bins,
            aes(x = bin_index, y = log2ratio + yoffset, label = name),
            color = "red", size = 2.8, fontface = "bold") +
  geom_segment(data = segments_df,
               aes(x = x_start, xend = x_end, y = seg_value, yend = seg_value),
               inherit.aes = FALSE, linewidth = 0.7, color = "darkorange") +
  geom_vline(data = chr_bounds %>% filter(chromosome != 1),
             aes(xintercept = min_bin - 0.5),
             linetype = "dashed", linewidth = 0.3, color = "grey60") +
  scale_x_continuous(name = "chromosome",
                     breaks = chr_labels$label_pos,
                     labels = chr_labels$chromosome,
                     expand = expansion(mult = c(0, 0))) +
  ylab("log2 ratio") +
  coord_cartesian(ylim = c(-3, 6)) +
  ggtitle(plot_title) +
  scale_color_manual(name = "CGHcall (bins)",
                     values = c(loss = "blue", gain = "red", amp = "red", neutral = "grey30")) +
  scale_alpha_continuous(range = c(0.2, 1), guide = "none") +
  annotate("text", x = x_min, y = Inf, label = left_label, hjust = -0.1, vjust = 1.2, size = 3) +
  annotate("text", x = x_max, y = Inf, label = right_label, hjust = 1.1, vjust = 1.2, size = 3) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

cnv_outfile <- file.path(output_dir, paste0(sample_id, "_CNV.png"))
safe_ggsave(cnv_outfile, p, width = 10, height = 4, dpi = 200)
message("Saved CNV plot: ", cnv_outfile)

segment_calls_outfile <- file.path(output_dir, paste0(sample_id, "_CGHcall_segments.tsv"))
segment_calls <- segments_df %>%
  filter(call_class != "neutral") %>%
  arrange(chrom, start) %>%
  transmute(sample = sample_id,
            chromosome = chrom,
            start = start,
            end = end,
            n_bins = n_bins,
            seg_mean_log2 = seg_value_raw,
            call = call_class,
            seg_call_code = seg_call,
            seg_prob = seg_prob)
write_tsv(segment_calls, segment_calls_outfile)
message("Saved CGHcall segment table: ", segment_calls_outfile)

qc_outfile <- file.path(output_dir, paste0(sample_id, "_QDNAseq_QC.tsv"))
qc_df <- data.frame(sample = sample_id,
                    expected_sigma = obs_sigma,
                    delta_sigma = delta_sigma,
                    called_segments = n_segments,
                    stringsAsFactors = FALSE)
write_tsv(qc_df, qc_outfile)

if (generate_locuszoom) {
  locus_dir <- file.path(output_dir, "locuszooms")
  dir.create(locus_dir, recursive = TRUE, showWarnings = FALSE)
  message("Generating locuszoom plots for ", sample_id)
  for (g in unique(genes_df$name)) {
    pz <- plot_gene_region_qc(cn, genes_df, g, sample_id)
    outfile <- file.path(locus_dir, paste0(sample_id, "_", g, ".png"))
    safe_ggsave(outfile, pz, width = 6, height = 4, dpi = 200)
  }
  message("Saved locuszooms to: ", locus_dir)
}

message("Single-sample CGHcall annotation completed successfully.")
