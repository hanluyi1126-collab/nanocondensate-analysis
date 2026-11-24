# Plot settings and exports shared by the two figure entry points.
configure_figure_theme <- function() {
  theme_set(theme_classic(base_size = 9, base_family = "Helvetica") +
    theme(axis.line = element_line(linewidth = 0.35),
          axis.ticks = element_line(linewidth = 0.3),
          axis.text = element_text(colour = "black"),
          plot.title = element_text(size = 10, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.key.width = grid::unit(9, "mm"),
          plot.margin = margin(3, 4, 3, 3, "mm"),
          strip.background = element_blank()))
}

save_figure_files <- function(plot, name, height_mm, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- function(ext) file.path(output_dir, paste0(name, ".", ext))
  tagged <- (plot + plot_annotation(tag_levels = "A")) &
    theme(plot.tag = element_text(face = "bold", size = 11))

  # Close each device even if patchwork fails during drawing.
  draw_and_close <- function() {
    on.exit(grDevices::dev.off())
    print(tagged)
  }

  svglite::svglite(path("svg"), width = 6.5, height = height_mm / 25.4)
  draw_and_close()
  grDevices::pdf(path("pdf"), width = 6.5, height = height_mm / 25.4,
                 family = "Helvetica", useDingbats = FALSE)
  draw_and_close()
  ragg::agg_png(path("png"), width = 165, height = height_mm, units = "mm", res = 300)
  draw_and_close()
  ragg::agg_tiff(path("tiff"), width = 165, height = height_mm, units = "mm",
                 res = 600, compression = "lzw")
  draw_and_close()
}

draw_droplet_figures <- function(data_dir, save_figure) {
  fusion <- read.csv(file.path(data_dir, "fusion_verified.csv"))
  potential_levels <- c(0, 10, 25)
  fusion$potential <- factor(fusion$potential_mV, levels = potential_levels,
                             labels = paste(potential_levels, "mV"))
  neutral <- subset(fusion, potential_mV == 0)
  colours <- c("#333333", "#347D9E", "#885080")

  encounters <- ggplot(neutral, aes(radius_um, encounter_s)) +
    geom_line(colour = colours[2], linewidth = 0.65) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "Droplet radius (µm)", y = expression(nu[enc]~(s^{-1})),
         title = "Brownian encounters")
  capillary <- ggplot(neutral, aes(radius_um, capillary_kBT)) +
    geom_line(colour = colours[2], linewidth = 0.65) + scale_x_log10() +
    labs(x = "Droplet radius (µm)", y = expression(Delta*G[gamma]/(k[B]*T)),
         title = "Capillary free energy")
  energy <- ggplot(fusion, aes(radius_um, barrier_kBT, colour = potential)) +
    geom_line(linewidth = 0.65) + scale_x_log10() + scale_colour_manual(values = colours) +
    labs(x = "Droplet radius (µm)", y = expression(U[el](0)/(k[B]*T)),
         title = "Formal HHF contact energy") +
    theme(legend.position = "top")
  contact <- ggplot(fusion, aes(radius_um, log10_contact_rate_s, colour = potential)) +
    geom_line(linewidth = 0.65) + scale_x_log10() + scale_colour_manual(values = colours) +
    labs(x = "Droplet radius (µm)", y = expression(log[10](nu[contact] / s^{-1})),
         title = "Formal contact arrival rate") +
    theme(legend.position = "top")
  save_figure((encounters + capillary) / (energy + contact), "Fig_S1", 125)

  state <- read.csv(file.path(data_dir, "state_grid.csv"))
  state$state <- factor(state$state,
    levels = c("Dispersed", "Coarsening droplets", "Nanocondensate window", "Gel/solid"),
    labels = c("Dispersed", "Coarsening", "Candidate window", "Network arrest"))
  state_map <- ggplot(state, aes(sticker_fraction, charged_fraction, fill = state)) +
    geom_raster() +
    scale_fill_manual(values = c("#E0E3E5", "#DFC18C", "#86B5AA", "#9B7C9F"), drop = FALSE) +
    coord_fixed() + labs(x = "Sticker fraction", y = "Charged fraction",
                         title = "Illustrative state map") +
    theme(legend.position = "bottom", legend.text = element_text(size = 7)) +
    guides(fill = guide_legend(nrow = 2))

  # Equations S11 and S12. The surface-to-volume ratio for diameter d is 6/d.
  diameter_nm <- exp(seq(log(50), log(2000), length.out = 120))
  density_per_nm3 <- 0.005
  charge_per_nm2 <- 0.03
  charges_per_molecule <- 5
  footprint_nm2 <- 10
  charge_fraction <- 6 * charge_per_nm2 / (charges_per_molecule * density_per_nm3 * diameter_nm)
  coverage_fraction <- 6 / (footprint_nm2 * density_per_nm3 * diameter_nm)
  budget <- rbind(data.frame(diameter_nm, fraction = charge_fraction, criterion = "Charge budget"),
                  data.frame(diameter_nm, fraction = coverage_fraction, criterion = "Complete coverage"))
  budgets <- ggplot(budget, aes(diameter_nm, fraction, colour = criterion)) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3) + geom_line(linewidth = 0.7) +
    scale_x_log10() + scale_y_log10() + scale_colour_manual(values = colours[2:3]) +
    labs(x = "Droplet diameter (nm)", y = "Required molecular fraction",
         title = "Interfacial molecule budgets") +
    theme(legend.position = "bottom", legend.text = element_text(size = 7))
  save_figure(state_map + budgets, "Fig_S2", 98)
}
