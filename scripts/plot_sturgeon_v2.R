#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 3){
  stop("Usage: Rscript plot_new_classifier.R <input_csv> <sample_name> <output_png>")
}

in_csv <- args[1]
sample_name <- args[2]
out_png <- args[3]

# ---------------- Load Data ----------------
df <- read.csv(in_csv, header = TRUE)

# Last column is 'probes'
probes_used <- df$probes[1]

# Reshape for plotting
df_plot <- df %>%
  select(-probes) %>%
  pivot_longer(cols = everything(), names_to = "tumor_class", values_to = "confidence")

# ---------------- Plot ----------------
p <- ggplot(df_plot, aes(x = tumor_class, y = confidence)) +
  geom_bar(stat = "identity", fill = "steelblue") +

  # Confidence threshold lines
  geom_hline(yintercept = 0.80, linetype = "dotted", color = "darkgoldenrod", linewidth = 0.7) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "red", linewidth = 0.7) +

  # Titles
  labs(
    title = sample_name,
    subtitle = paste("Number of probes:", probes_used),
    x = "Tumor Class",
    y = "Confidence Score"
  ) +

  # Force y-axis to start at 0
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +

  # Theme adjustments
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), # vertical x-axis labels
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    plot.title = element_text(face = "bold")
  )

# ---------------- Save PNG ----------------
ggsave(filename = out_png, plot = p, width = 10, height = 5, dpi = 300)
