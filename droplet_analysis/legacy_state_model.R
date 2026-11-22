#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

out_dir <- "outputs/state"
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
grid <- subset(grid, sticker_fraction + charged_fraction <= 0.85)

# 1) Condensation: non-electrostatic stickers dominate; residual opposite-charge
# pairing contributes weakly to cohesion.
grid$cohesion_index <- with(grid,
  4.0 * sticker_fraction +
  1.2 * charged_fraction * (1 - charge_asymmetry^2) * screening
)

# 2) Anti-fusion: net surface charge, modular segregation and ionic screening.
grid$surface_charge_index <- with(grid,
  5.0 * charged_fraction * charge_asymmetry * charge_blockiness * screening
)
grid$barrier_kBT <- 200 * grid$surface_charge_index^2

# 3) Solidification: sticker connectivity plus internal complementary-charge
# crosslinking. High charge asymmetry reduces complementary ionic pairing.
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
  geom_contour_filled(
    aes(z = state_code),
    breaks = c(0.5, 1.5, 2.5, 3.5, 4.5),
    show.legend = TRUE
  ) +
  geom_contour(aes(z = cohesion_index), breaks = 1, color = "#3C4650",
               linewidth = 0.35, linetype = "dashed") +
  geom_contour(aes(z = barrier_kBT), breaks = c(10, 20), color = "white",
               linewidth = 0.35) +
  geom_contour(aes(z = network_index), breaks = 1, color = "black",
               linewidth = 0.45) +
  scale_fill_manual(
    values = unname(cols),
    labels = names(cols),
    drop = FALSE
  ) +
  coord_fixed(xlim = c(0, 0.6), ylim = c(0, 0.6), expand = FALSE) +
  labs(
    title = "Predicted nanocondensate formation window",
    subtitle = "Dashed: condensation threshold; white: 10/20 kBT fusion barriers; black: gelation threshold",
    x = "Non-electrostatic sticker fraction",
    y = "Charged-component fraction",
    fill = "Predicted state"
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

cat(sprintf("Predicted target window: sticker %.3f--%.3f; charge %.3f--%.3f\n",
            min(target$sticker_fraction), max(target$sticker_fraction),
            min(target$charged_fraction), max(target$charged_fraction)))
