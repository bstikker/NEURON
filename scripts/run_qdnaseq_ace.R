# run_qdnaseq_ace.R
# Usage: Rscript run_qdnaseq_ace.R <sample_id> <adjusted_bam_path>

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

if (length(args) < 2) {
  stop("Not enough arguments provided. Usage: Rscript run_qdnaseq_ace.R <sample_id> <adjusted_bam_path>")
}

sample_id <- args[1]
bam_path <- normalizePath(args[2], mustWork = TRUE)
cat(sprintf("Running QDNAseq + ACE for %s\n", sample_id))

base_dir <- file.path(repo_root, "results", sample_id, "QDNAseq_ACE")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

suppressMessages({
  library(QDNAseq)
  library(DNAcopy)
  library(ACE)
  library(ggplot2)
  library(Biobase)
})

bins_path <- file.path(repo_root, "reference", "bins", "T2T.500kbp.SR50.rda")
if (!file.exists(bins_path)) stop("QDNAseq bins file not found: ", bins_path)
load(bins_path)
if (!exists("bins")) stop("Object 'bins' not found after loading ", bins_path)

readCounts <- binReadCounts(bins, bamfiles = bam_path)
readCountsFiltered <- applyFilters(estimateCorrection(readCounts))
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersSmooth <- smoothOutlierBins(normalizeBins(copyNumbers))
copyNumbersSegmented <- normalizeSegmentedBins(segmentBins(copyNumbersSmooth, transformFun = "sqrt"))

rds_path <- file.path(base_dir, paste0(sample_id, "_copyNumbersSegmented.rds"))
saveRDS(copyNumbersSegmented, file = rds_path)

bed_path <- file.path(base_dir, paste0(sample_id, "_500kbp.bed"))
write.table(
  as.data.frame(fData(copyNumbersSegmented)),
  file = bed_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

seg_path <- file.path(base_dir, paste0(sample_id, "_500kbp.seg"))
write.table(
  as.data.frame(Biobase::assayDataElement(copyNumbersSegmented, "segmented")),
  file = seg_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

seg_png <- file.path(base_dir, paste0(sample_id, "_segmented.png"))
png(seg_png, width = 1200, height = 400)
plot(copyNumbersSegmented)
dev.off()

summary_file <- file.path(base_dir, "ACE_summary.tsv")
plot_file <- file.path(base_dir, paste0(sample_id, "_ACE_matrixplot.png"))

ace_fun_candidates <- c(
  file.path(repo_root, "scripts", "ACE_functions.R"),
  file.path(repo_root, "scripts", "R_scripts", "ACE_functions.R")
)
ace_fun_path <- ace_fun_candidates[file.exists(ace_fun_candidates)][1]
if (!is.na(ace_fun_path) && nzchar(ace_fun_path)) {
  source(ace_fun_path)
}

if (!exists("squaremodel")) {
  stop("squaremodel() was not found. Add ACE_functions.R to scripts/ or scripts/R_scripts/.")
}

ACE_result <- squaremodel(copyNumbersSegmented, QDNAseqobjectsample = TRUE, penalty = 0.5, penploidy = 0.5)
best_fit <- ACE_result$errordf[order(ACE_result$errordf$error), ][1, ]
write.table(best_fit, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
ggsave(filename = plot_file, plot = ACE_result$matrixplot, width = 8, height = 6, dpi = 300)

cat("QDNAseq + ACE analysis completed successfully.\n")
