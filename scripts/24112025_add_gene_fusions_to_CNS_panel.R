#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

in_bed  <- file.path(repo_root, "reference", "genes", "relevant_genes_with_chm13v2_500kb_bin_nrs.bed")
out_bed <- file.path(repo_root, "reference", "genes", "relevant_genes_with_chm13v2_500kb_bin_nrs_fusions_singlebin.bed")

genes <- read_tsv(
  in_bed,
  col_names = c("chrom", "start", "end", "name", "binnr", "yoffset"),
  col_types = cols(
    chrom   = col_character(),
    start   = col_integer(),
    end     = col_integer(),
    name    = col_character(),
    binnr   = col_integer(),
    yoffset = col_double()
  ),
  trim_ws = TRUE
)

fusion_map <- c(
  "CDKN2A" = "CDKN2A/B",
  "CDKN2B" = "CDKN2A/B",
  "FGFR3"  = "FGFR3/TACC3",
  "TACC3"  = "FGFR3/TACC3",
  "FGFR1"  = "FGFR1-TACC1",
  "TACC1"  = "FGFR1-TACC1"
)

genes2 <- genes %>%
  dplyr::mutate(fusion_name = dplyr::recode(name, !!!fusion_map, .default = name))

genes3 <- genes2 %>%
  dplyr::mutate(overlap_bp = end - start)

fused_single_bin <- genes3 %>%
  dplyr::group_by(fusion_name) %>%
  dplyr::slice_max(overlap_bp, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(chrom, start) %>%
  dplyr::mutate(name = fusion_name) %>%
  dplyr::select(chrom, start, end, name, binnr, yoffset)

write_tsv(fused_single_bin, out_bed, col_names = TRUE)
cat("Wrote fusion-aware single-bin BED with", nrow(fused_single_bin), "rows to:\n", out_bed, "\n")
