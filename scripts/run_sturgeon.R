# run_sturgeon.R (enhanced)
# Usage: Rscript run_sturgeon.R <bam_path> <sample_id>

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript run_sturgeon.R <bam_path> <sample_id>")
}

bam_file <- normalizePath(args[1], mustWork = TRUE)
sample_id <- args[2]

cat(sprintf("Running Sturgeon pipeline for sample: %s\n", sample_id))

output_dir <- file.path(repo_root, "results", sample_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

adjusted_bam <- file.path(output_dir, "adjusted_merged_sorted.bam")
modkit_out <- file.path(output_dir, "modkit_extracted.txt")
sturgeon_shell <- file.path(repo_root, "scripts", "sturgeon_guppy_modkit_SLURM.sh")

# ---------------------------
# Step 0: Estimate RAM based on BAM size
# ---------------------------
bam_size_gb <- file.info(bam_file)$size / 1024^3
estimated_mem_gb <- max(20, min(60, ceiling(bam_size_gb * 5)))  # simple heuristic: 5GB per GB BAM
cat(sprintf("BAM size: %.2f GB, recommended RAM: ~%d GB\n", bam_size_gb, estimated_mem_gb))

# ---------------------------
# Step 1: Adjust mods
# ---------------------------
cmd_adjust <- sprintf("modkit adjust-mods --convert h m %s %s", shQuote(bam_file), shQuote(adjusted_bam))
cat(sprintf("Running: %s\n", cmd_adjust))
status <- system(cmd_adjust)
if (status != 0) stop("modkit adjust-mods failed.")

# ---------------------------
# Step 1b: Index BAM if missing (use samtools)
# ---------------------------
index_file <- paste0(adjusted_bam, ".bai")
if (!file.exists(index_file)) {
  cat("Index file missing. Creating BAM index with samtools...\n")
  cmd_index <- sprintf("samtools index %s", shQuote(adjusted_bam))
  status <- system(cmd_index)
  if (status != 0) stop("samtools index failed.")
} else {
  cat("Found BAM index.\n")
}

# ---------------------------
# Step 2: Extract modkit data
# ---------------------------
threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset="4"))
cat(sprintf("Using %d threads for extraction\n", threads))

cmd_extract <- sprintf("modkit extract full -t %d %s %s", threads, shQuote(adjusted_bam), shQuote(modkit_out))
cat(sprintf("Running: %s\n", cmd_extract))
status <- system(cmd_extract)
if (status != 0) stop("modkit extract failed.")

# ---------------------------
# Step 3: Run Sturgeon prediction
# ---------------------------
cmd_sturgeon <- sprintf("%s %s", shQuote(sturgeon_shell), shQuote(output_dir))
cat(sprintf("Running: %s\n", cmd_sturgeon))
status <- system(cmd_sturgeon)
if (status != 0) stop("Sturgeon prediction shell script failed.")

cat("Sturgeon pipeline completed successfully.\n")
