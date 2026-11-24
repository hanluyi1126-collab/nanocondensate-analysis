#!/usr/bin/env Rscript
# Check Eqs. S1-S10 against the archived figure inputs in Supplementary Data 1.
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!file.exists(script_path)) script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
module_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
data_dir <- file.path(module_dir, "data", "drop")
fusion <- read.csv(file.path(data_dir, "fusion_verified.csv"))

# Reconstruct the 180 log-spaced radii and three potentials in Table S1.
radii_um <- exp(seq(log(0.05), log(2), length.out = 180))
samples <- expand.grid(potential_mV = c(0, 10, 25), radius_um = radii_um)
fusion <- fusion[order(fusion$radius_um, fusion$potential_mV), ]
stopifnot(nrow(fusion) == nrow(samples),
          all(fusion$potential_mV == samples$potential_mV),
          max(abs(fusion$radius_um - samples$radius_um)) < 1e-14)

thermal_energy <- 1.380649e-23 * 298.15
viscosity <- 0.890e-3
permittivity <- 78.3 * 8.8541878128e-12
kappa <- sqrt(2 * (1.602176634e-19)^2 * 6.02214076e23 * 1000 * 1e-7 /
              (permittivity * thermal_energy))
radius_m <- samples$radius_um * 1e-6
diffusion <- thermal_energy / (6 * pi * viscosity * radius_m)
number_density <- 0.01 / ((4 / 3) * pi * radius_m^3)
encounter <- 16 * pi * diffusion * radius_m * number_density
capillary <- 4 * pi * 1e-6 * radius_m^2 * (2^(2/3) - 2) / thermal_energy
amplitude <- 2 * pi * permittivity * radius_m * (samples$potential_mV / 1000)^2 /
  thermal_energy
contact_energy <- amplitude * log(2)

# Eq. S8 uses y = 2R/r. Subtract the contact energy before exponentiation
# so even the 917 kBT endpoint can be integrated in double precision.
midpoints <- (seq_len(24000) - 0.5) / 24000
integrate_pair <- function(i) {
  if (samples$potential_mV[i] == 0) return(c(adaptive = 0, midpoint = 0))
  scaled_integrand <- function(y) {
    gap <- 2 * radius_m[i] * (1 - y) / y
    energy <- amplitude[i] * log1p(exp(-kappa * gap))
    exp(energy - contact_energy[i])
  }
  adaptive <- integrate(scaled_integrand, 0, 1, abs.tol = 1e-12, rel.tol = 1e-10)$value
  midpoint <- mean(scaled_integrand(midpoints))
  (contact_energy[i] + log(c(adaptive = adaptive, midpoint = midpoint))) / log(10)
}
log_w <- vapply(seq_len(nrow(samples)), integrate_pair, c(adaptive = 0, midpoint = 0))
rate <- log10(encounter) - log_w["adaptive", ]
midpoint_difference <- log_w["midpoint", ] - log_w["adaptive", ]

check_values <- function(label, observed, expected, tolerance) {
  err <- max(abs(observed - expected))
  if (!is.finite(err) || err > tolerance) {
    stop(sprintf("%s: maximum absolute difference %.5g exceeds %.5g", label, err, tolerance))
  }
  message(sprintf("%-32s max. absolute difference %.4g", label, err))
}
check_values("diameter (nm)", fusion$diameter_nm, radius_m * 2e9, 1e-10)
check_values("diffusivity (um^2/s)", fusion$diffusion_um2_s, diffusion * 1e12, 1e-12)
check_values("encounter rate (1/s)", fusion$encounter_s, encounter, 1e-11)
check_values("capillary energy (kBT)", fusion$capillary_kBT, capillary, 1e-10)
check_values("HHF contact energy (kBT)", fusion$barrier_kBT, contact_energy, 1e-11)
check_values("kappa times radius", fusion$kappa_R, kappa * radius_m, 1e-12)
check_values("adaptive log10 W", fusion$log10_W, log_w["adaptive", ], 1e-10)
check_values("adaptive log10 contact rate", fusion$log10_contact_rate_s, rate, 1e-10)
check_values("midpoint minus adaptive log10 W", fusion$log10_midpoint_minus_adaptive,
             midpoint_difference, 1e-10)
stopifnot(all(rate <= log10(encounter) + 1e-12))
message(sprintf("Largest integration-method difference: %.7g in log10 W", max(abs(midpoint_difference))))

grid <- read.csv(file.path(data_dir, "state_grid.csv"))
s <- grid$sticker_fraction
c_frac <- grid$charged_fraction
positions <- data.frame(s = round(s * 200), c = round(c_frac * 200))
stopifnot(!anyDuplicated(positions), nrow(grid) == 12106,
          all(abs(s * 200 - positions$s) < 1e-10), all(abs(c_frac * 200 - positions$c) < 1e-10),
          all(positions >= 0), all(positions <= 120), all(positions$s + positions$c <= 170))
complete_grid <- expand.grid(s = 0:120, c = 0:120)
complete_grid <- subset(complete_grid, s + c <= 170)
membership <- merge(complete_grid, transform(positions, archived = TRUE), all.x = TRUE)
missing <- subset(membership, is.na(archived))
stopifnot(nrow(missing) == 50, all(missing$s + missing$c == 170))

cohesion <- 4 * s + 1.2 * c_frac * (1 - 0.6^2) * 0.5
charge <- 5 * c_frac * 0.6 * 0.8 * 0.5
barrier <- 200 * charge^2
network <- 3 * s^2 + 6 * c_frac^2 * (1 - 0.6^2)
check_values("cohesion score", grid$cohesion_index, cohesion, 1e-12)
check_values("charge score", grid$surface_charge_index, charge, 1e-12)
check_values("repulsion score", grid$barrier_kBT, barrier, 1e-12)
check_values("network score", grid$network_index, network, 1e-12)

classify <- function(multiplier) {
  ifelse(cohesion < multiplier, 1L,
    ifelse(network >= multiplier, 4L, ifelse(barrier < 10 * multiplier, 2L, 3L)))
}
state_names <- c("Dispersed", "Coarsening droplets", "Nanocondensate window", "Gel/solid")
nominal <- classify(1)
stopifnot(all(grid$state_code == nominal), all(grid$state == state_names[nominal]))
sensitivity <- vapply(c(0.8, 1, 1.2), function(t) sum(classify(t) == 3L), integer(1))
stopifnot(all(sensitivity == c(2100, 2318, 2486)))
message(sprintf("Table S3 matches all %d archived classifications. The inclusive grid has %d additional boundary points.",
                nrow(grid), nrow(missing)))
