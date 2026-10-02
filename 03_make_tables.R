# Manuscript: Does a mortality schedule need a Makeham term? Calibrating the likelihood-ratio test
# Read saved simulations; create the manuscript table and supporting summaries.
# This script performs no fitting and draws no new random numbers.
script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project_dir <- if (length(script_argument)) dirname(normalizePath(sub("^--file=", "", script_argument[1]))) else getwd()
source(file.path(project_dir, "01_functions.R"))
folder <- results_folder(project_dir)
manifest <- readRDS(file.path(folder, "manifest.rds"))
results <- readRDS(file.path(folder, "replicate_results.rds"))
limit <- readRDS(file.path(folder, "joint_limit.rds"))
settings <- manifest$signature$settings
# Read the mode from the saved run, not the current console environment.
mode <- if (settings$B == 10000L) "full" else if (settings$B == 20L) "smoke" else "pilot"
info <- manifest$information
scenarios <- manifest$scenarios
out <- file.path(folder, "tables")
dir.create(out, showWarnings = FALSE)
cut_boundary <- qchisq(1 - 2 * settings$alpha, 1)
cut_usual <- qchisq(1 - settings$alpha, 1)

# --- "Finite-sample calibration" and "Local power and detectability" --------
# Every probability includes the number of usable fits and its Monte Carlo SE.
# If a fit failed, the reported estimate is conditional on successful fits.
# failure_lower/upper show what the rate could be across ALL replicates.
# Split once instead of scanning the entire results table for every scenario.
# Column sums summarize several test rules together; no per-replicate data frames.
# Function: summarize_events
# Each column is a test outcome or zero-mass indicator for the same scenario.
# Summarize all columns together, retaining Monte Carlo uncertainty and the
# range of probabilities possible if failed fits were resolved either way.
summarize_events <- function(events) {
  total <- nrow(events)
  n <- colSums(!is.na(events))
  yes <- colSums(events, na.rm = TRUE)
  p <- ifelse(n > 0, yes / n, NA_real_)
  z <- qnorm(.975)
  centre <- (p + z^2/(2*n))/(1 + z^2/n)
  half <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/(1 + z^2/n)
  data.frame(metric = colnames(events), total = total, successful = n,
    failed = total - n, estimate = p, mcse = sqrt(p*(1-p)/n),
    mc_lower = centre - half, mc_upper = centre + half,
    failure_lower = yes/total, failure_upper = (yes + total - n)/total,
    row.names = NULL)
}
# Group the rows once; reuse those groups for all test rules.
by_scenario <- split(seq_len(nrow(results)), results$scenario)
summary <- do.call(rbind, lapply(scenarios$scenario, function(s) {
  d <- results[by_scenario[[as.character(s)]], ]
  T <- d$GM_T
  gm <- summarize_events(cbind(mass_zero = T == 0,
                               boundary = T > cut_boundary, usual = T > cut_usual))
  gm$model <- "GM"
  if (d$delta[1] == 0) {
    T <- d$joint_T
    joint <- summarize_events(cbind(mass_zero = T == 0, boundary = T > cut_boundary,
                                    usual = T > cut_usual, joint_critical = T > limit$critical))
    joint$model <- "joint"
    gm <- rbind(gm, joint)
  }
  cbind(scenarios[rep(match(s, scenarios$scenario), nrow(gm)), ],
        gm[, c("model", "metric", "total", "successful", "failed", "estimate", "mcse",
               "mc_lower", "mc_upper", "failure_lower", "failure_upper")])
}))
rownames(summary) <- NULL
write.csv(summary, file.path(out, "all_probability_summaries.csv"), row.names = FALSE)
saveRDS(summary, file.path(out, "all_probability_summaries.rds"))
calibration <- subset(summary, model == "GM" & delta == 0)
power <- subset(summary, model == "GM" & metric %in% c("boundary", "usual"))
threshold <- subset(power, is_threshold)
joint <- subset(summary, model == "joint")
write.csv(calibration, file.path(out, "table1_calibration_details.csv"), row.names = FALSE)
write.csv(power, file.path(out, "power_curves.csv"), row.names = FALSE)
write.csv(threshold, file.path(out, "detection_threshold.csv"), row.names = FALSE)
write.csv(joint, file.path(out, "joint_boundary_summary.csv"), row.names = FALSE)

# The limit itself is simulated: report its Monte Carlo uncertainty too.
limit_summary <- do.call(rbind, lapply(c("mass_zero", "boundary", "usual", "joint_critical"),
  function(metric) {
    event <- switch(metric, mass_zero = limit$T0 == 0,
      boundary = limit$T0 > cut_boundary, usual = limit$T0 > cut_usual,
      joint_critical = limit$T0 > limit$critical)
    cbind(metric = metric, probability_summary(event))
  }))
# The last row uses the same draws that determined the quantile; it is NOT an
# independent check of calibration. Finite-sample simulations are independent.
write.csv(limit_summary, file.path(out, "joint_limit_summary.csv"), row.names = FALSE)
# Collect the information and cutoffs used in the manuscript text.
constants <- data.frame(j_eff = info$j, k11 = info$K[1, 1], k12 = info$K[1, 2],
  k22 = info$K[2, 2], rho = info$rho, exact_joint_mass_zero = info$mass_zero,
  boundary_critical = cut_boundary, usual_critical = cut_usual,
  simulated_joint_critical = limit$critical, gaussian_draws = length(limit$T0),
  threshold_delta = qnorm(1 - settings$alpha) + qnorm(settings$target_power))
write.csv(constants, file.path(out, "numbers_for_manuscript.csv"), row.names = FALSE)

# --- One LaTeX table for the manuscript --------------------------------------
# Failures stay in the table. Use the CSV for MC SEs, Wilson intervals and bounds.
# This is a standalone snippet, not an edit to makeham_research_notes.tex.
# Function: value
# Retrieve and format one calibration-table entry. Missing entries appear as
# "--" so an unavailable result cannot look like a zero probability.
value <- function(n, metric, field = "estimate") {
  z <- calibration[calibration$n == n & calibration$metric == metric, field]
  if (!length(z) || is.na(z)) return("--")
  if (field == "failed") as.character(z) else sprintf("%.3f", z)
}
tol_parts <- strsplit(format(settings$controls$zero_tolerance, scientific = TRUE), "e", fixed = TRUE)[[1]]
tol_tex <- paste0(tol_parts[1], "\\times10^{", as.integer(tol_parts[2]), "}")
lines <- c("% Generated from saved Monte Carlo results by 03_make_tables.R.",
  paste0("% Mode: ", mode, "."),
  "\\begin{table}[ht]", "\\centering",
  "\\caption{Finite-sample calibration of the Makeham likelihood-ratio test.}",
  "\\label{tab:sim_calibration}", "\\begin{tabular}{rrrrr}", "\\hline",
  "$n$ & $\\Pr(T\\approx0)$ & $\\Pr(T>2.7055)$ & $\\Pr(T>3.8415)$ & Failed \\\\",
  "\\hline")
# Table 1 writes the four manuscript exposure totals as powers of ten.
for (n in settings$exposure_totals) {
  lines <- c(lines, paste(paste0("$10^{", as.integer(log10(n)), "}$"),
    value(n, "mass_zero"), value(n, "boundary"), value(n, "usual"),
    value(n, "mass_zero", "failed"), sep = " & "))
  lines[length(lines)] <- paste0(lines[length(lines)], " \\\\")
}
lines <- c(lines, "\\hline", "Asymptotic & 0.500 & 0.050 & 0.025 & -- \\\\",
  "\\hline", "\\end{tabular}", "\\par\\smallskip",
  paste0("{\\footnotesize $B=", format(settings$B, big.mark = "{,}", scientific = FALSE, trim = TRUE),
    "$ per exposure level. Probabilities use successful fits; failures are shown separately. ",
    "We classify $|T_{\\mathrm{raw}}|\\leq ", tol_tex,
    "$ as numerical zero. ", if (mode != "full") "Pilot/check run; not final manuscript results." else "", "}"), "\\end{table}")
# Caption thresholds above match alpha=.05. Stop rather than silently mislabel.
if (settings$alpha != .05) stop("Update the LaTeX headings for your chosen alpha.")
writeLines(lines, file.path(out, "table1_calibration.tex"))

# --- Parameter estimates and optimization diagnostics -----------------------
# Checkpoints retain all starts, fitted age-specific rates and raw death counts.
# This convenient CSV extracts the selected fit for each model and replicate.
# Function: extract_parameters
# Read one checkpoint and collect each available selected fit into a compact
# data frame. Return parameter estimates, fit quality and optimizer diagnostics.
# Only one detailed checkpoint is read at a time; the compact rows are combined
# and written once after extraction.
extract_parameters <- function(job) {
  saved <- readRDS(file.path(folder, "checkpoints", sprintf("chunk_%05d.rds", job$id)))
  entries <- unlist(lapply(seq_along(saved$records), function(i) {
    fits <- saved$records[[i]]$fits
    models <- intersect(c("G", "GM", "GG", "GGM"), names(fits))
    models <- models[vapply(fits[models], function(fit) !is.null(fit$parameters), logical(1))]
    lapply(models, function(model)
      list(replicate = job$first + i - 1L, model = model, fit = fits[[model]]))
  }), recursive = FALSE)
  if (!length(entries)) return(NULL)
  parameters <- t(vapply(entries, function(z) unname(z$fit$parameters), numeric(4)))
  colnames(parameters) <- c("lambda", "beta", "kappa", "gamma")
  data.frame(scenario = job$scenario,
    replicate = vapply(entries, function(z) z$replicate, numeric(1)),
    model = vapply(entries, function(z) z$model, character(1)),
    n = saved$scenario$n, delta = saved$scenario$delta, parameters,
    ok = vapply(entries, function(z) z$fit$ok, logical(1)),
    log_likelihood = vapply(entries, function(z) z$fit$log_likelihood, numeric(1)),
    half_deviance = vapply(entries, function(z) z$fit$objective, numeric(1)),
    projected_score = vapply(entries, function(z) z$fit$score, numeric(1)),
    selected_start = vapply(entries, function(z) z$fit$source, character(1)),
    attempts = vapply(entries, function(z) nrow(z$fit$attempts), integer(1)),
    message = vapply(entries, function(z) z$fit$message, character(1)))
}
parameter_tables <- lapply(manifest$jobs, extract_parameters)
parameters <- do.call(rbind, parameter_tables)
if (!is.null(parameters)) {
  rownames(parameters) <- NULL
  write.csv(parameters, file.path(out, "fitted_parameters.csv"), row.names = FALSE)
}
rm(parameter_tables, parameters)
# Assess sensitivity of the estimated zero mass to a purely numerical choice.
zero_checks <- do.call(rbind, lapply(settings$exposure_totals, function(n) {
  d <- results[results$delta == 0 & results$n == n, ]
  do.call(rbind, lapply(c("GM", "joint"), function(model) {
    raw <- d[[paste0(model, "_raw_T")]]
    do.call(rbind, lapply(c(1e-9, 1e-7, 1e-5), function(tolerance)
      cbind(n = n, model = model, tolerance = tolerance,
            probability_summary(abs(raw) <= tolerance))))
  }))
}))
write.csv(zero_checks, file.path(out, "zero_tolerance_sensitivity.csv"), row.names = FALSE)
if (any(summary$failed > 0)) warning("Some fits failed. Review failed counts, bounds and optimizer diagnostics before interpreting rates.")
print(constants, row.names = FALSE)
print(calibration[, c("n", "metric", "estimate", "mcse", "failed")], row.names = FALSE)
message("Tables and supporting CSVs saved in: ", out)
