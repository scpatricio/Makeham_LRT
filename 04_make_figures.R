# Make the two manuscript figures from SAVED simulation results using ggplot2.
# Run 03_make_tables.R first. No fitting or simulation takes place here.
# Saves two vector PDFs in results/figures/. Includes legends; no caption files.
# If needed, install once: install.packages("ggplot2")
library(ggplot2)

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
default_project_dir <- if (length(script_argument)) dirname(normalizePath(sub("^--file=", "", script_argument[1]))) else getwd()
project_dir <- Sys.getenv("MAKEHAM_PROJECT_DIR", unset = default_project_dir)
source(file.path(project_dir, "01_functions.R"))
folder <- results_folder(project_dir)
manifest <- readRDS(file.path(folder, "manifest.rds"))
settings <- manifest$signature$settings
summary <- readRDS(file.path(folder, "tables", "all_probability_summaries.rds"))
results <- readRDS(file.path(folder, "replicate_results.rds"))
limit <- readRDS(file.path(folder, "joint_limit.rds"))
out <- Sys.getenv("MAKEHAM_FIGURE_DIR", unset = file.path(folder, "figures"))
dir.create(out, showWarnings = FALSE)

# Use the same exposure colours in both figures, in exposure_totals order:
# blue, red, gold and green for the four manuscript exposure totals.
exposures <- as.character(settings$exposure_totals)
colours <- setNames(rep(c("#0077B6", "#E63946", "#FFB703", "#2A9D8F"),
                        length.out = length(exposures)), exposures)
shapes <- setNames(rep(c(15, 16, 17, 18), length.out = length(exposures)), exposures)
# Compact mathematical labels for the manuscript legend: N = 10^4, ..., N = 10^7.
exposure_exponents <- round(log10(settings$exposure_totals))
exposure_labels <- parse(text = paste0("N == 10^", exposure_exponents))
cut_boundary <- qchisq(1 - 2 * settings$alpha, 1)
threshold_delta <- qnorm(1 - settings$alpha) + qnorm(settings$target_power)
pilot_note <- if (settings$B < 10000L) "Pilot run - not final manuscript results" else NULL


# Horizontal legend adjustment used in the manuscript figures.
h_shift = 15

figure_theme <- theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.key.width = grid::unit(1, "cm"),
    legend.key.height = grid::unit(0.2, "cm"),
    legend.spacing.y = grid::unit(1, "pt"),
    legend.box.spacing = grid::unit(2, "pt"),
    legend.margin = margin(0, h_shift, 0, -h_shift),
    
    plot.title.position = "plot",
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0
    ),
    
    plot.margin = margin(0, 0, 0, 0)
  )


# --- "Local power and detectability" ---------------------------------------
# Points: simulated rejection probabilities under the corrected cutoff.
# Optional error bars: uncomment geom_errorbar below to show pointwise 95%
# Wilson Monte Carlo intervals. They are off in the supplied manuscript style.
# Black solid/dashed curves: boundary/conventional local-power approximations.
# The exact threshold point stays in the saved tables; it is omitted here
# because it nearly coincides with delta=2.5 in the regular plotting grid.
power_data <- summary[summary$model == "GM" & summary$metric == "boundary" &
                        !summary$is_threshold, ]
power_data$exposure <- factor(as.character(power_data$n), levels = exposures)
power_grid <- data.frame(delta = seq(0, max(settings$delta_grid), length.out = 401))
power_grid$boundary <- pnorm(power_grid$delta - qnorm(1 - settings$alpha))
power_grid$usual <- pnorm(power_grid$delta - qnorm(1 - settings$alpha/2))

boundary_approx_label <- "Boundary approx."
usual_approx_label <- "Conventional approx."
target_label <- paste0(format(100 * settings$target_power, trim = TRUE), "% target")
figure_power <- ggplot() +
  geom_hline(aes(yintercept = settings$target_power, linetype = target_label), colour = "grey65", key_glyph = "path", linewidth = .8) +
  geom_vline(xintercept = threshold_delta, colour = "grey65", linetype = "dotted", show.legend = FALSE, linewidth = .8) +
  geom_line(data = power_grid, aes(delta, boundary, linetype = boundary_approx_label), colour = "black", linewidth = .9) +
  geom_line(data = power_grid, aes(delta, usual, linetype = usual_approx_label), colour = "black",
            linewidth = .9) +
  # geom_errorbar(data = power_data,
  #               aes(delta, ymin = mc_lower, ymax = mc_upper, colour = exposure),
  #               width = .05, alpha = .5, linewidth = .45, show.legend = FALSE) +
  geom_line(data = power_data, aes(delta, estimate, colour = exposure, group = exposure),
            linewidth = .6) +
  geom_point(data = power_data, aes(delta, estimate, colour = exposure, shape = exposure),
             size = 2.3) +
  scale_colour_manual(name = "", values = colours, breaks = exposures, labels = exposure_labels) +
  scale_shape_manual(name = "", values = shapes, breaks = exposures, labels = exposure_labels) +
  scale_linetype_manual(name = "", values = setNames(c("solid", "dashed", "dotted"),
                                                     c(boundary_approx_label, usual_approx_label, target_label)),
                        breaks = c(boundary_approx_label, usual_approx_label, target_label)) +
  guides(colour = guide_legend(order = 1, ncol = 2, byrow = TRUE),
         shape = guide_legend(order = 1, ncol = 2, byrow = TRUE),
         linetype = guide_legend(order = 2, ncol = 2, byrow = TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, .2)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Local power of the Makeham likelihood-ratio test",
       subtitle = pilot_note,
       x = expression(delta == u * sqrt(j[eff])), y = "Rejection probability") +
  figure_theme

# --- "The joint boundary at zero frailty" ----------------------------------
# Upper-tail probabilities preserve the mass at zero; no density smoothing.
# At threshold t, each curve gives the probability of a statistic above t.
tail_grid <- seq(0, max(8, unname(quantile(limit$T0, .995))), length.out = 501)
limit_data <- data.frame(t = tail_grid,
                         joint = 1 - ecdf(limit$T0)(tail_grid),
                         half_chisq = .5 * pchisq(tail_grid, 1, lower.tail = FALSE))
# Each exposure gets its empirical tail from successful null fits only.
tail_data <- do.call(rbind, lapply(settings$exposure_totals, function(n) {
  T <- results$joint_T[results$delta == 0 & results$n == n]
  T <- T[is.finite(T)]
  if (!length(T)) return(NULL)
  data.frame(t = tail_grid, probability = 1 - ecdf(T)(tail_grid),
             exposure = factor(as.character(n), levels = exposures))
}))

figure_joint <- ggplot() +
  geom_vline(aes(xintercept = cut_boundary, linetype = "Boundary cutoff"), colour = "grey55", key_glyph = "path", linewidth = .8) +
  geom_vline(aes(xintercept = limit$critical, linetype = "Joint cutoff"), colour = "grey55", key_glyph = "path", linewidth = .8) +
  geom_line(data = limit_data, aes(t, joint, linetype = "Joint limit"), colour = "black", linewidth = 1) +
  geom_line(data = limit_data, aes(t, half_chisq, linetype = "Half-chi-square ref."), colour = "black",
            linewidth = .9) +
  scale_colour_manual(name = "", values = colours, breaks = exposures, labels = exposure_labels) +
  scale_linetype_manual(name = "", values = c("Joint limit" = "solid",
                                              "Half-chi-square ref." = "dashed",
                                              "Boundary cutoff" = "dotted",
                                              "Joint cutoff" = "dotdash"),
                        breaks = c("Joint limit", "Half-chi-square ref.",
                                   "Boundary cutoff", "Joint cutoff")) +
  guides(colour = guide_legend(order = 1, ncol = 2, byrow = TRUE),
         linetype = guide_legend(order = 2, ncol = 2, byrow = TRUE)) +
  coord_cartesian(ylim = c(0, max(.65, limit_data$joint[1] + .02))) +
  labs(title = "Makeham test with zero true frailty variance",
       subtitle = pilot_note,
       x = "Likelihood-ratio threshold t", y = expression(Pr(T > t))) +
  figure_theme
if (!is.null(tail_data)) {
  figure_joint <- figure_joint +
    geom_line(data = tail_data, aes(t, probability, colour = exposure, group = exposure),
              linewidth = .65)
}

# --- Save PDFs only ---------------------------------------------------------
# The ggplot objects remain available to edit or display in RStudio.
ggsave(file.path(out, "figure1_local_power.pdf"), plot = figure_power,
       width = 5, height = 4, units = "in", bg = "white")
ggsave(file.path(out, "figure2_joint_boundary.pdf"), plot = figure_joint,
       width = 5, height = 4, units = "in", bg = "white")
if (interactive()) {
  print(figure_power)
  print(figure_joint)
}
message("Saved the two ggplot2 PDFs to: ", out)

