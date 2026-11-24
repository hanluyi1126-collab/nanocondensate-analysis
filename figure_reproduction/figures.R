#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg)
# Rscript 4.4 can encode spaces in its --file argument.
if (!file.exists(script_path)) script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
source(file.path(script_dir, "..", "droplet_analysis", "plot_helpers.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript figures.R [figure_reproduction_directory]")
root <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else script_dir
membrane_dir <- file.path(root, "data", "membrane")
output_dir <- file.path(root, "figures")

# Several membrane tables have empty columns that are unused in these panels.
# Check only the measurements contributing to a panel before passing them to ggplot.
read_series <- function(name, n, numeric_columns, other_columns = character()) {
  path <- file.path(membrane_dir, paste0(name, ".csv"))
  dat <- read.csv(path, check.names = FALSE)
  missing <- setdiff(c(numeric_columns, other_columns), names(dat))
  if (length(missing)) stop(name, ": missing column(s) ", paste(missing, collapse = ", "))
  if (nrow(dat) != n) stop(name, ": expected ", n, " archived rows, found ", nrow(dat))
  for (col in numeric_columns) {
    values <- dat[[col]]
    if (!is.numeric(values) || any(!is.finite(values))) {
      stop(name, ": expected finite numeric values in ", col)
    }
  }
  dat
}

diameter <- read_series("diameter", 201, c("diameter_nm", "centre_tension_pN_per_nm",
  "maximum_tension_pN_per_nm", "maximum_tension_position_fraction", "force_residual"))
if (!setequal(diameter$diameter_nm, 50:250)) {
  stop("diameter: expected one solution per diameter from 50 to 250 nm")
}
refinement <- read_series("convergence", 35, c("diameter_nm", "membrane_elements",
  "centre_tension_pN_per_nm", "force_residual"))
mesh_sizes <- c(240, 480, 960, 1920, 3840)
mesh_diameters <- c(60, 80, 94, 100, 120, 160, 200)
mesh_counts <- with(refinement, table(factor(diameter_nm, levels = mesh_diameters),
                                      factor(membrane_elements, levels = mesh_sizes)))
if (any(mesh_counts != 1L)) stop("convergence: incomplete or duplicated refinement series")
finest <- subset(refinement, membrane_elements == max(mesh_sizes))
coarser <- subset(refinement, membrane_elements < max(mesh_sizes))
reference <- finest$centre_tension_pN_per_nm[match(coarser$diameter_nm, finest$diameter_nm)]
if (any(reference == 0)) stop("convergence: zero reference tension in relative error")
coarser$centre_error_pct <- 100 * abs(coarser$centre_tension_pN_per_nm / reference - 1)

grid <- read_series("cohesive_grid", 1554, "force_residual")
ensemble <- read_series("ensemble", 2048, "force_residual")
residuals <- data.frame(force_residual = c(diameter$force_residual, refinement$force_residual,
                                           grid$force_residual, ensemble$force_residual))
if (any(residuals$force_residual <= 0) || any(coarser$centre_error_pct <= 0)) {
  stop("Fig. S3: the archived log axes require positive residuals and relative errors")
}

trajectories <- read_series("lipid_replicates", 24, c("radius_nm", "force_over_2pi_pN"))
if (any(table(factor(trajectories$radius_nm, levels = c(2, 3, 4, 5, 6, 8))) != 4L)) {
  stop("lipid_replicates: expected four trajectory means per archived radius")
}
fit <- lm(force_over_2pi_pN ~ radius_nm, data = trajectories, na.action = na.fail)
bootstrap <- read_series("bootstrap_verified", 20000, "line_tension_pN")
line_tension_interval <- quantile(bootstrap$line_tension_pN, c(0.025, 0.975), type = 7)

phases <- read_series("phase_boundaries", 336, character(),
                      c("classification", "first_50_percent_crossing_nm"))
outcome_labels <- c(
  "pore at both diameter limits" = "Pore at both limits",
  "small-particle wrap; large-particle pore" = "Small wrap / large pore",
  "small-particle pore; large-particle wrap" = "Small pore / large wrap",
  "wrap at both diameter limits" = "Wrap at both limits"
)
if (any(!phases$classification %in% names(outcome_labels))) {
  stop("phase_boundaries: missing or unrecognised endpoint classification")
}
counts <- as.data.frame(table(factor(phases$classification, levels = names(outcome_labels))))
names(counts) <- c("classification", "count")
counts$label <- unname(outcome_labels[as.character(counts$classification)])
directional <- subset(phases, classification == "small-particle pore; large-particle wrap")
if (!nrow(directional) || !is.numeric(directional$first_50_percent_crossing_nm) ||
    any(!is.finite(directional$first_50_percent_crossing_nm))) {
  stop("phase_boundaries: expected a finite first crossing for each small-pore / large-wrap branch")
}

curves <- read_series("selected_probability_curves", 48, c("diameter_nm",
  "two_stage_pore_probability_before_wrap", "lysis_reference_time_ms"))
curve_diameters <- c(seq(50, 110, 5), 120, 130, 150)
reference_times <- c(0.001, 0.01, 0.1)
curve_counts <- with(curves, table(factor(diameter_nm, levels = curve_diameters),
                                   factor(lysis_reference_time_ms, levels = reference_times)))
if (any(curve_counts != 1L)) stop("selected_probability_curves: incomplete or duplicated curves")
probability <- curves$two_stage_pore_probability_before_wrap
if (any(probability < 0 | probability > 1)) stop("selected_probability_curves: invalid probability")
curves$reference_time <- factor(curves$lysis_reference_time_ms, levels = reference_times)

configure_figure_theme()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Cannot create figure directory: ", output_dir)
save_figure <- function(plot, name, height_mm) save_figure_files(plot, name, height_mm, output_dir)
draw_droplet_figures(file.path(root, "data", "drop"), save_figure)

blue <- "#347D9E"
purple <- "#885080"
stress <- rbind(
  data.frame(diameter_nm = diameter$diameter_nm, tension = diameter$centre_tension_pN_per_nm,
             site = "Centre"),
  data.frame(diameter_nm = diameter$diameter_nm, tension = diameter$maximum_tension_pN_per_nm,
             site = "Spatial maximum")
)
stress_plot <- ggplot(stress, aes(diameter_nm, tension, colour = site)) +
  geom_line(linewidth = 0.65) +
  geom_hline(yintercept = 8, linetype = 2, linewidth = 0.3) +
  scale_colour_manual(values = c(blue, purple)) +
  labs(x = "Diameter (nm)", y = expression(tau~(pN~nm^{-1})), title = "Centre and peak tension") +
  theme(legend.position = "top")
peak_position <- ggplot(diameter, aes(diameter_nm, maximum_tension_position_fraction)) +
  geom_line(colour = purple, linewidth = 0.65) +
  expand_limits(y = c(0, 1.15)) +
  labs(x = "Diameter (nm)", y = expression(r[max]/a[c]), title = "Position of peak tension")
mesh_plot <- ggplot(coarser, aes(membrane_elements, centre_error_pct, group = diameter_nm)) +
  geom_line(colour = blue, alpha = 0.6, linewidth = 0.45) +
  geom_point(size = 1.2, colour = blue) +
  scale_x_log10(breaks = head(mesh_sizes, -1)) + scale_y_log10() +
  labs(x = "Membrane elements", y = "Centre error versus 3,840 (%)", title = "Mesh refinement")
residual_plot <- ggplot(residuals, aes(force_residual)) +
  stat_ecdf(geom = "step", colour = blue, linewidth = 0.6) +
  geom_vline(xintercept = 1e-6, linetype = 2, linewidth = 0.3) + scale_x_log10() +
  labs(x = "Force residual (pN)", y = "Cumulative fraction",
       title = sprintf("Force residuals (n = %s)", format(nrow(residuals), big.mark = ",")))
save_figure((stress_plot + peak_position) / (mesh_plot + residual_plot), "Fig_S3", 125)

force_plot <- ggplot(trajectories, aes(radius_nm, force_over_2pi_pN)) +
  geom_point(colour = blue, alpha = 0.7, size = 1.8) +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], colour = purple, linewidth = 0.6) +
  labs(x = "Restrained pore radius (nm)", y = expression(F[hold]/(2*pi)~(pN)),
       title = sprintf("Trajectory means (n = %d)", nrow(trajectories)))
interval_plot <- ggplot(bootstrap, aes(line_tension_pN)) +
  geom_histogram(bins = 45, fill = blue, colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = line_tension_interval, linetype = 2, linewidth = 0.35) +
  labs(x = expression(lambda~(pN)), y = "Bootstrap resamples", title = "Line-tension uncertainty")
save_figure(force_plot + interval_plot, "Fig_S4", 87)

endpoint_plot <- ggplot(counts, aes(reorder(label, count), count)) +
  geom_col(fill = blue, width = 0.65) + geom_text(aes(label = count), hjust = -0.2, size = 3) +
  coord_flip(clip = "off") + scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(x = NULL, y = "Parameter combinations",
       title = "Endpoint classes") +
  theme(axis.text.y = element_text(size = 7))
crossing_plot <- ggplot(directional, aes(first_50_percent_crossing_nm)) +
  geom_histogram(binwidth = 6, fill = purple, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = median(directional$first_50_percent_crossing_nm), linetype = 2,
             linewidth = 0.4) +
  labs(x = "First 50% crossing diameter (nm)", y = "Directional branches",
       title = sprintf("Small pore / large wrap (n = %d)", nrow(directional)))
probability_plot <- ggplot(curves, aes(diameter_nm, two_stage_pore_probability_before_wrap,
                                       colour = reference_time)) +
  geom_line(linewidth = 0.65) + geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3) +
  scale_colour_manual(name = "Reference time (ms)", values = c(blue, purple, "#9A773C")) +
  labs(x = "Diameter (nm)", y = "Conditional pore probability", title = "Defect reference time") +
  theme(legend.position = "right", legend.title = element_text(size = 8))
save_figure((endpoint_plot + crossing_plot) / probability_plot, "Fig_S5", 130)

writeLines(capture.output(sessionInfo()), file.path(root, "R_session_info.txt"))
