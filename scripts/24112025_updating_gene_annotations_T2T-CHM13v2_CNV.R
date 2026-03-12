#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

gtf_file <- file.path(repo_root, "reference", "GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf")
out_bed  <- file.path(repo_root, "reference", "CNS_gene_panel.bed")

genes_raw <- c(
  "MDM4", "PTEN", "MGMT", "CCND1", "CCND2", "CDK4", "MDM2", "RB1", "TP53",
  "NF1", "PPM1D", "MYCN", "GLI2", "FGFR3/TACC3", "PDGFRA", "TERT", "MYB",
  "EGFR", "CDK6", "MET", "BRAF", "FGFR1-TACC1", "MYBL1", "MYC", "CDKN2A/B",
  "PTCH1"
)

fusion_map <- list(
  "FGFR3/TACC3" = c("FGFR3", "TACC3"),
  "FGFR1-TACC1" = c("FGFR1", "TACC1"),
  "CDKN2A/B"    = c("CDKN2A", "CDKN2B")
)

expanded_genes <- unlist(
  lapply(genes_raw, function(g) if (g %in% names(fusion_map)) fusion_map[[g]] else g),
  use.names = FALSE
) %>% unique()

gtf <- read_tsv(
  gtf_file,
  comment = "#",
  col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
  col_types = cols(
    seqname   = col_character(),
    source    = col_character(),
    feature   = col_character(),
    start     = col_integer(),
    end       = col_integer(),
    score     = col_character(),
    strand    = col_character(),
    frame     = col_character(),
    attribute = col_character()
  )
)

genes_gtf <- gtf %>% filter(feature == "gene")

extract_gene_symbol <- function(attr) {
  m1 <- str_match(attr, 'gene_name "([^"]+)"')
  sym <- m1[, 2]
  need_gene <- is.na(sym)
  if (any(need_gene)) {
    m2 <- str_match(attr[need_gene], 'gene "([^"]+)"')
    sym[need_gene] <- m2[, 2]
  }
  sym
}

genes_gtf <- genes_gtf %>%
  mutate(gene_symbol = extract_gene_symbol(attribute)) %>%
  filter(!is.na(gene_symbol))

genes_sel <- genes_gtf %>%
  filter(gene_symbol %in% expanded_genes)

write_tsv(
  genes_sel %>% select(seqname, start, end, gene_symbol, strand),
  out_bed,
  col_names = TRUE
)

cat("Wrote CNS gene panel BED to:\n", out_bed, "\n")
