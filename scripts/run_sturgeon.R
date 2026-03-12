# run_sturgeon.R
# Usage: Rscript run_sturgeon.R <bam_path> <sample_id>

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

if (length(args) < 2) {
  stop("Not enough arguments provided. Usage: Rscript run_sturgeon.R <bam_path> <sample_id>")
}

bam_file <- normalizePath(args[1], mustWork = TRUE)
sample_id <- args[2]

cat(sprintf("Running Sturgeon pipeline for sample: %s\n", sample_id))

output_dir <- file.path(repo_root, "results", sample_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

adjusted_bam <- file.path(output_dir, "adjusted_merged_sorted.bam")
modkit_out <- file.path(output_dir, "modkit_extracted.txt")
sturgeon_modkit_shellscript <- file.path(repo_root, "scripts", "sturgeon_guppy_modkit_SLURM.sh")

cmd_adjust <- sprintf(
  "modkit adjust-mods --convert h m %s %s",
  shQuote(bam_file), shQuote(adjusted_bam)
)
cat(sprintf("Running: %s\n", cmd_adjust))
status <- system(cmd_adjust)
if (status != 0) stop("modkit adjust-mods failed.")

cmd_extract <- sprintf(
  "modkit extract full -t 10 %s %s",
  shQuote(adjusted_bam), shQuote(modkit_out)
)
cat(sprintf("Running: %s\n", cmd_extract))
status <- system(cmd_extract)
if (status != 0) stop("modkit extract failed.")

cmd_sturgeon <- sprintf(
  "%s %s",
  shQuote(sturgeon_modkit_shellscript),
  shQuote(output_dir)
)
cat(sprintf("Running: %s\n", cmd_sturgeon))
status <- system(cmd_sturgeon)
if (status != 0) stop("Sturgeon prediction shell script failed.")

cat("Sturgeon pipeline completed successfully.\n")
