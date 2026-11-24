suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(svglite)
  library(ragg)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!file.exists(script_path)) script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
module_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
out_dir <- file.path(module_dir, "outputs", "fusion")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Supplementary Table S1. SI units throughout, except where the unit is in the name.
kB <- 1.380649e-23
T <- 298.15
kBT <- kB * T
eta <- 0.890e-3
eps0 <- 8.8541878128e-12
epsr <- 78.3
eps <- eps0 * epsr
N_A <- 6.02214076e23
e <- 1.602176634e-19
ionic_strength_M <- 1e-7
kappa <- sqrt(2 * e^2 * N_A * 1000 * ionic_strength_M / (eps * kBT))
lambda_D <- 1 / kappa
phi <- 0.01
gamma <- 1e-6
K_brown <- 8 * kBT / (3 * eta)

radius_um <- exp(seq(log(0.05), log(2), length.out = 180))
radius_m <- radius_um * 1e-6

base_df <- data.frame(
  radius_um = radius_um,
  D_um2_s = kBT / (6 * pi * eta * radius_m) * 1e12,
  number_density_um3 = 3 * phi / (4 * pi * radius_m^3) / 1e18,
  encounter_s = K_brown * (3 * phi / (4 * pi * radius_m^3)),
  deltaG_cap_kBT = 4 * pi * gamma * radius_m^2 * (2^(2/3) - 2) / kBT
)

barrier_kBT <- function(R_m, psi_V) {
  2 * pi * eps * R_m * psi_V^2 * log(2) / kBT
}

# Fuchs stability ratio using the HHF constant-potential double-layer potential.
# Substitution y = 2R/r maps the semi-infinite integral to y in (0, 1).
log_fuchs_W <- function(R_m, psi_V, n_grid = 24000) {
  if (psi_V == 0) return(0)
  y <- (seq_len(n_grid) - 0.5) / n_grid
  x <- 2 * (1 - y) / y
  pref <- 2 * pi * eps * R_m * psi_V^2 / kBT
  u <- pref * log1p(exp(-kappa * R_m * x))
  m <- max(u)
  m + log(mean(exp(u - m)))
}

psi_levels <- c(0, 0.010, 0.025)
psi_labels <- c("Neutral", "10 mV", "25 mV")
fusion_list <- vector("list", length(psi_levels))
for (j in seq_along(psi_levels)) {
  psi <- psi_levels[j]
  lw <- vapply(radius_m, log_fuchs_W, numeric(1), psi_V = psi)
  fusion_list[[j]] <- data.frame(
    radius_um = radius_um,
    psi = factor(psi_labels[j], levels = psi_labels),
    barrier_kBT = barrier_kBT(radius_m, psi),
    log10_fusion_s = log10(base_df$encounter_s) - lw / log(10)
  )
}
fusion_df <- do.call(rbind, fusion_list)
# These historical column names assume unit post-contact fusion probability.
# The display floor applies only to this diagnostic plot.
fusion_df$log10_fusion_display <- pmax(fusion_df$log10_fusion_s, -50)

write.csv(
  merge(base_df, fusion_df, by = "radius_um", all = TRUE),
  file.path(out_dir, "pure_water_droplet_fusion_source_data.csv"),
  row.names = FALSE
)

COL <- c(
  teal = "#2A8C8C",
  teal_light = "#BFE3DF",
  blue = "#185AA9",
  blue_light = "#D7E7F6",
  coral = "#D85B45",
  grey = "#6C737A",
  grey_light = "#E8EBED",
  ink = "#202326"
)

theme_set(
  theme_classic(base_size = 7.2, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = COL["ink"]),
      axis.ticks = element_line(linewidth = 0.35, colour = COL["ink"]),
      axis.title = element_text(size = 7.3),
      axis.text = element_text(size = 6.6, colour = COL["ink"]),
      legend.title = element_blank(),
      legend.text = element_text(size = 6.5),
      legend.key.width = unit(4.5, "mm"),
      legend.position = "top",
      plot.title = element_text(size = 7.5, face = "bold", margin = margin(b = 2)),
      plot.margin = margin(3, 5, 3, 4),
      panel.grid = element_blank()
    )
)

circle_df <- function(cx, cy, r, id, n = 120) {
  th <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(th), y = cy + r * sin(th), id = id)
}

droplets <- do.call(rbind, list(
  circle_df(0.95, 1.55, 0.22, "s1"),
  circle_df(1.65, 2.12, 0.22, "s2"),
  circle_df(2.05, 1.25, 0.22, "s3"),
  circle_df(2.85, 1.78, 0.22, "s4"),
  circle_df(4.30, 1.62, 0.36, "c1"),
  circle_df(5.22, 1.62, 0.36, "c2"),
  circle_df(7.50, 1.62, 0.72, "l1"),
  circle_df(9.15, 1.62, 0.72, "l2")
))

charge_pts <- do.call(rbind, lapply(list(c(4.30, 1.62), c(5.22, 1.62)), function(z) {
  th <- seq(0, 2 * pi, length.out = 9)[-9]
  data.frame(x = z[1] + 0.45 * cos(th), y = z[2] + 0.45 * sin(th))
}))

p_a <- ggplot() +
  annotate("rect", xmin = 0.25, xmax = 3.35, ymin = 0.45, ymax = 2.72,
           fill = "#F4FAF9", colour = NA) +
  annotate("rect", xmin = 3.65, xmax = 5.86, ymin = 0.45, ymax = 2.72,
           fill = "#F3F7FC", colour = NA) +
  annotate("rect", xmin = 6.25, xmax = 9.95, ymin = 0.45, ymax = 2.72,
           fill = "#F7F8F9", colour = NA) +
  geom_polygon(data = droplets, aes(x, y, group = id),
               fill = COL["teal_light"], colour = COL["teal"], linewidth = 0.55) +
  geom_point(data = charge_pts, aes(x, y), shape = 21, size = 2.2,
             stroke = 0.45, fill = "white", colour = COL["blue"]) +
  geom_text(data = charge_pts, aes(x, y, label = "+"),
            size = 2.1, family = "Helvetica", colour = COL["blue"]) +
  annotate("segment", x = 4.00, xend = 3.58, y = 1.62, yend = 1.62,
           linewidth = 0.75, colour = COL["blue"],
           arrow = arrow(length = unit(1.7, "mm"), type = "closed")) +
  annotate("segment", x = 5.52, xend = 5.94, y = 1.62, yend = 1.62,
           linewidth = 0.75, colour = COL["blue"],
           arrow = arrow(length = unit(1.7, "mm"), type = "closed")) +
  annotate("curve", x = 0.62, y = 1.13, xend = 0.90, yend = 0.80,
           curvature = 0.45, linetype = 2, linewidth = 0.35, colour = COL["grey"],
           arrow = arrow(length = unit(1.2, "mm"))) +
  annotate("curve", x = 1.47, y = 2.38, xend = 1.12, yend = 2.55,
           curvature = -0.45, linetype = 2, linewidth = 0.35, colour = COL["grey"],
           arrow = arrow(length = unit(1.2, "mm"))) +
  annotate("curve", x = 2.28, y = 1.03, xend = 2.62, yend = 0.76,
           curvature = 0.5, linetype = 2, linewidth = 0.35, colour = COL["grey"],
           arrow = arrow(length = unit(1.2, "mm"))) +
  annotate("segment", x = 7.02, xend = 6.82, y = 2.35, yend = 2.48,
           linetype = 2, linewidth = 0.35, colour = COL["grey"],
           arrow = arrow(length = unit(1.1, "mm"))) +
  annotate("segment", x = 9.62, xend = 9.80, y = 0.94, yend = 0.82,
           linetype = 2, linewidth = 0.35, colour = COL["grey"],
           arrow = arrow(length = unit(1.1, "mm"))) +
  annotate("text", x = 1.80, y = 2.57, label = "Small radii (50-200 nm)",
           family = "Helvetica", fontface = "bold", size = 2.8, colour = COL["ink"]) +
  annotate("text", x = 1.80, y = 0.61,
           label = "Fast diffusion, frequent encounters",
           family = "Helvetica", size = 2.35, colour = COL["grey"]) +
  annotate("text", x = 4.76, y = 2.53, label = "Charged interface",
           family = "Helvetica", fontface = "bold", size = 2.8, colour = COL["ink"]) +
  annotate("text", x = 4.76, y = 0.61,
           label = "Repulsive approach potential",
           family = "Helvetica", size = 2.35, colour = COL["blue"]) +
  annotate("text", x = 8.33, y = 2.57, label = "Large radii (0.5-2 µm)",
           family = "Helvetica", fontface = "bold", size = 2.8, colour = COL["ink"]) +
  annotate("text", x = 8.33, y = 0.61,
           label = "Slow diffusion, rare encounters",
           family = "Helvetica", size = 2.35, colour = COL["grey"]) +
  annotate("segment", x = 0.55, xend = 9.55, y = 2.96, yend = 2.96,
           linewidth = 0.55, colour = COL["ink"],
           arrow = arrow(length = unit(1.7, "mm"), type = "closed")) +
  annotate("text", x = 5.05, y = 3.15,
           label = "Droplet radius increases from 50 nm to 2 µm",
           family = "Helvetica", size = 2.6, colour = COL["ink"]) +
  coord_cartesian(xlim = c(0.2, 10), ylim = c(0.35, 3.30), clip = "off") +
  labs(title = "a   Droplet encounters in the idealized water model") +
  theme_void(base_family = "Helvetica", base_size = 7.2) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 2)),
    plot.margin = margin(5, 4, 0, 4)
  )

size_breaks <- c(0.05, 0.1, 0.2, 0.5, 1, 2)

p_b <- ggplot(base_df, aes(radius_um, encounter_s)) +
  geom_line(linewidth = 0.85, colour = COL["teal"]) +
  geom_point(data = subset(base_df, radius_um %in% radius_um[round(seq(1, nrow(base_df), length.out = 6))]),
             size = 1.4, colour = COL["teal"]) +
  scale_x_log10(breaks = size_breaks, labels = size_breaks) +
  scale_y_log10() +
  annotate("text", x = 0.11, y = 0.7, label = expression(nu[enc] %prop% R^{-3}),
           family = "Helvetica", size = 2.5, colour = COL["teal"]) +
  labs(title = "b   Brownian encounters", x = "Droplet radius, R (µm)",
       y = expression("Encounter frequency, "*nu[enc]*" (s"^{-1}*")"))

bar_df <- subset(fusion_df, psi != "Neutral")
p_c <- ggplot(bar_df, aes(radius_um, barrier_kBT, colour = psi)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = COL["grey"]) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("10 mV" = unname(COL["teal"]), "25 mV" = unname(COL["blue"]))) +
  scale_x_log10(breaks = size_breaks, labels = size_breaks) +
  annotate("text", x = 0.065, y = 875, label = "Barrier > 0 for all R",
           family = "Helvetica", size = 2.35, hjust = 0, colour = COL["coral"]) +
  labs(title = "c   Formal HHF contact energy", x = "Droplet radius, R (µm)",
       y = expression("Positive barrier, "*Delta*G[el]^ddagger*" / "*k[B]*T)) +
  theme(legend.position = "top")

p_d <- ggplot(base_df, aes(radius_um, deltaG_cap_kBT)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = COL["grey"]) +
  geom_line(linewidth = 0.85, colour = COL["coral"]) +
  scale_x_log10(breaks = size_breaks, labels = size_breaks) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(-5000, -1000, -100, -10, -1, 0)
  ) +
  annotate("text", x = 0.07, y = -1300, label = expression(Delta*G[gamma] < 0),
           family = "Helvetica", size = 2.5, hjust = 0, colour = COL["coral"]) +
  labs(title = "d   Capillary free-energy change", x = "Droplet radius, R (µm)",
       y = expression(Delta*G[gamma]*" / "*k[B]*T))

p_e <- ggplot(fusion_df, aes(radius_um, log10_fusion_display, colour = psi)) +
  annotate("rect", xmin = 0.05, xmax = 2, ymin = -50, ymax = -20,
           fill = COL["grey_light"], alpha = 0.55) +
  geom_hline(yintercept = -20, linetype = 2, linewidth = 0.35, colour = COL["grey"]) +
  geom_line(linewidth = 0.85) +
  scale_colour_manual(values = c(
    "Neutral" = unname(COL["ink"]),
    "10 mV" = unname(COL["teal"]),
    "25 mV" = unname(COL["blue"])
  )) +
  scale_x_log10(breaks = size_breaks, labels = size_breaks) +
  scale_y_continuous(breaks = c(-50, -40, -20, 0, 2), limits = c(-50, 2.8)) +
  annotate("text", x = 0.055, y = -35, label = "Low model rates",
           family = "Helvetica", size = 2.35, hjust = 0, colour = COL["grey"]) +
  annotate("text", x = 1.85, y = -48, label = "Display floor at -50",
           family = "Helvetica", size = 2.15, hjust = 1, colour = COL["grey"]) +
  labs(title = "e   Idealized contact arrival", x = "Droplet radius, R (µm)",
       y = expression(log[10](nu[contact]/s^{-1}))) +
  theme(legend.position = "top")

bottom_row <- p_b + p_c + p_d + plot_layout(widths = c(1, 1, 1))
final_plot <- p_a / bottom_row / p_e +
  plot_layout(heights = c(1.15, 1.0, 1.0)) +
  plot_annotation(
    caption = paste0(
      "Assumptions: 25 °C, idealized water (Debye length = ", sprintf("%.2f", lambda_D * 1e6),
      " µm), volume fraction 1%, interfacial tension 1 µN/m. Constant-potential HHF/Fuchs model; no van der Waals, steric, hydration or film-drainage terms."
    ),
    theme = theme(
      plot.caption = element_text(family = "Helvetica", size = 5.8, colour = COL["grey"], hjust = 0),
      plot.margin = margin(3, 4, 3, 4)
    )
  )

base <- file.path(out_dir, "pure_water_droplet_fusion_nature")
width_mm <- 183
height_mm <- 145
w <- width_mm / 25.4
h <- height_mm / 25.4

svglite::svglite(paste0(base, ".svg"), width = w, height = h)
print(final_plot)
dev.off()

grDevices::pdf(file = paste0(base, ".pdf"), width = w, height = h,
               family = "Helvetica", useDingbats = FALSE)
print(final_plot)
dev.off()

ragg::agg_png(paste0(base, ".png"), width = width_mm, height = height_mm,
              units = "mm", res = 600, background = "white")
print(final_plot)
dev.off()

ragg::agg_tiff(paste0(base, ".tiff"), width = width_mm, height = height_mm,
               units = "mm", res = 600, compression = "lzw", background = "white")
print(final_plot)
dev.off()

message("Saved: ", base)
