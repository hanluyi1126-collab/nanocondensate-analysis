#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!file.exists(script_path)) script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
module_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
out_dir <- file.path(module_dir, "outputs", "state")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Generic, dimensionless screening model. These are provisional assumptions,
# not fitted constants for a particular protein/polymer system.
charge_asymmetry <- 0.60  # |Q+ - Q-| / (Q+ + Q-)
charge_blockiness <- 0.80 # separation of cohesive and charged modules
salt_mM <- 150
screening <- 1 / sqrt(1 + salt_mM / 50)

grid <- expand.grid(
  sticker_fraction = seq(0, 0.60, by = 0.005),
  charged_fraction = seq(0, 0.60, by = 0.005)
)
# Use integer grid indices for the inclusive boundary in Supplementary Eq. S10.
# Retain the original floating-point selection to reproduce the archived figure.
grid$archived <- with(grid, sticker_fraction + charged_fraction <= 0.85)
grid <- subset(grid, round(200 * sticker_fraction) + round(200 * charged_fraction) <= 170)

# Prescribed cohesion, repulsion and connectivity scores (Supplementary Eq. S10).
grid$cohesion_index <- with(grid,
  4.0 * sticker_fraction +
  1.2 * charged_fraction * (1 - charge_asymmetry^2) * screening
)

grid$surface_charge_index <- with(grid,
  5.0 * charged_fraction * charge_asymmetry * charge_blockiness * screening
)
grid$barrier_kBT <- 200 * grid$surface_charge_index^2

grid$network_index <- with(grid,
  3.0 * sticker_fraction^2 +
  6.0 * charged_fraction^2 * (1 - charge_asymmetry^2)
)

grid$state <- with(grid, ifelse(
  cohesion_index < 1.0, "Dispersed",
  ifelse(network_index >= 1.0, "Gel/solid",
    ifelse(barrier_kBT < 10, "Coarsening droplets", "Nanocondensate window")
  )
))
grid$state <- factor(grid$state,
  levels = c("Dispersed", "Coarsening droplets", "Nanocondensate window", "Gel/solid"))
grid$state_code <- as.numeric(grid$state)

write.csv(subset(grid, !archived, select = -archived),
          file.path(out_dir, "additional_boundary_points.csv"), row.names = FALSE)
class_counts <- data.frame(state = levels(grid$state),
                           archived = as.integer(table(grid$state[grid$archived])),
                           inclusive_grid = as.integer(table(grid$state)))
write.csv(class_counts, file.path(out_dir, "grid_class_counts.csv"), row.names = FALSE)
grid <- subset(grid, archived, select = -archived)
write.csv(grid, file.path(out_dir, "nanocondensate_window_grid.csv"), row.names = FALSE)

target <- subset(grid, state == "Nanocondensate window")
summary_df <- data.frame(
  quantity = c("sticker_fraction_min", "sticker_fraction_max",
               "charged_fraction_min", "charged_fraction_max",
               "barrier_kBT_min", "barrier_kBT_median", "barrier_kBT_max"),
  value = c(min(target$sticker_fraction), max(target$sticker_fraction),
            min(target$charged_fraction), max(target$charged_fraction),
            min(target$barrier_kBT), median(target$barrier_kBT), max(target$barrier_kBT))
)
write.csv(summary_df, file.path(out_dir, "nanocondensate_window_summary.csv"), row.names = FALSE)

cols <- c("Dispersed" = "#DDE4E8", "Coarsening droplets" = "#F0BE70",
          "Nanocondensate window" = "#4C9B83", "Gel/solid" = "#76537D")

p <- ggplot(grid, aes(sticker_fraction, charged_fraction)) +
  geom_raster(aes(fill = state)) +
  geom_contour(aes(z = cohesion_index), breaks = 1, color = "#3C4650",
               linewidth = 0.35, linetype = "dashed") +
  geom_contour(aes(z = barrier_kBT), breaks = c(10, 20), color = "white",
               linewidth = 0.35) +
  geom_contour(aes(z = network_index), breaks = 1, color = "black",
               linewidth = 0.45) +
  scale_fill_manual(
    values = cols,
    labels = c("Dispersed", "Coarsening", "Candidate window", "Network arrest"),
    drop = FALSE
  ) +
  coord_fixed(xlim = c(0, 0.6), ylim = c(0, 0.6), expand = FALSE) +
  labs(
    title = "Illustrative composition map",
    subtitle = "Dashed: cohesion = 1; white: repulsion = 10, 20; black: network = 1",
    x = "Non-electrostatic sticker fraction",
    y = "Charged-component fraction",
    fill = "Assigned class"
  ) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.35),
        axis.ticks = element_line(linewidth = 0.35),
        plot.title = element_text(size = 9, face = "bold"),
        legend.title = element_text(size = 7), legend.text = element_text(size = 7),
        legend.position = "right")

ggsave(file.path(out_dir, "nanocondensate_window.png"), p,
       width = 145, height = 105, units = "mm", dpi = 400, device = ragg::agg_png)
grDevices::pdf(file.path(out_dir, "nanocondensate_window.pdf"),
               width = 145 / 25.4, height = 105 / 25.4,
               family = "Helvetica", useDingbats = FALSE,
               compress = TRUE, version = "1.4")
print(p)
dev.off()

cat(sprintf("Candidate grid range: sticker %.3f--%.3f; charge %.3f--%.3f\n",
            min(target$sticker_fraction), max(target$sticker_fraction),
            min(target$charged_fraction), max(target$charged_fraction)))
