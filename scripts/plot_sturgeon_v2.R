#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3){
  stop("Usage: Rscript plot_sturgeon_v2.R <sturgeon_v2_csv> <sample_id> <out_png>")
}

csv_file <- args[1]
sample_id <- args[2]
out_png <- args[3]

# Read the CSV
df <- read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE)

# Extract predicted classes and confidence
class_cols <- setdiff(colnames(df), c("probes","number_probes","sample","Sample_ID"))
scores <- as.numeric(df[1, class_cols])

plot_df <- data.frame(
  predicted_class = class_cols,
  display_class = dplyr::case_when(
    class_cols %in% c("A_IDH_HG","Glioma IDH - A IDH - HG") ~ "Glioma IDH - A IDH - HG",
    class_cols %in% c("O_IDH","Glioma IDH - O IDH - O IDH") ~ "Glioma IDH - O IDH - O IDH",
    class_cols %in% c("GBM_RTK1","Glioblastoma - GBM - RTK I") ~ "Glioblastoma - GBM - RTK I",
    class_cols %in% c("GBM_RTK2","Glioblastoma - GBM - RTK II") ~ "Glioblastoma - GBM - RTK II",
    class_cols %in% c("DMG_K27","Glioblastoma - DMG - K27") ~ "Glioblastoma - DMG - K27",
    TRUE ~ class_cols
  ),
  confidence = scores,
  stringsAsFactors = FALSE
)

# Shared palette for consistent colors
all_classes <- sort(unique(plot_df$display_class))
palette_vec <- setNames(hcl.colors(length(all_classes),"Dark 3"), all_classes)

# Plot
p <- ggplot(plot_df, aes(x=display_class, y=confidence, fill=display_class)) +
  geom_col() +
  scale_fill_manual(values = palette_vec) +
  ylim(0,1) +
  labs(title = paste0("Sturgeon v2 - ", sample_id),
       x = "Tumor class",
       y = "Confidence") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
        legend.position = "none",
        plot.title = element_text(face="bold"))

# Save PNG at smaller size to fit final_report
ggsave(filename = out_png, plot = p, width = 7, height = 4, dpi = 300)
cat("Sturgeon v2 plot written to", out_png, "\n")
