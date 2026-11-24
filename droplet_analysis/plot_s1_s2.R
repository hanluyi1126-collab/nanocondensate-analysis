#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!file.exists(script_path)) script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
module_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
source(file.path(module_dir, "plot_helpers.R"))
configure_figure_theme()
output_dir <- file.path(module_dir, "outputs", "figures")
draw_droplet_figures(file.path(module_dir, "data", "drop"), function(plot, name, height_mm) {
  save_figure_files(plot, name, height_mm, output_dir)
})
