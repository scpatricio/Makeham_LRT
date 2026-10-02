# Manuscript: Does a mortality schedule need a Makeham term? Calibrating the likelihood-ratio test
# This file contains the calculations. Run 02_run_simulations.R to use them.
# Base R only: no packages need to be installed.

# --- Common design: "Finite-sample calibration" ------------------------------
# Function: make_design
# Build one row per age with the null Gompertz rate and exposure proportion.
# The proportions sum to one; multiplying them by n gives person-years at each
# age. Survivorship is the manuscript design. The exponential option permits
# a separate sensitivity run with a user-specified decay rate.
make_design <- function(ages, lambda, beta, exposure_pattern = "survivorship",
                        exposure_decay = NULL) {
  stopifnot(length(ages) >= 4, all(diff(ages) > 0), min(ages) >= 0,
            lambda > 0, beta > 0)
  h <- lambda * exp(beta * ages)
  if (exposure_pattern == "survivorship") {
    log_weight <- -lambda / beta * (exp(beta * ages) - exp(beta * min(ages)))
  } else if (exposure_pattern == "exponential") {
    if (is.null(exposure_decay) || exposure_decay < 0)
      stop("Set a nonnegative exposure_decay for exponential exposures.")
    log_weight <- -exposure_decay * (ages - min(ages))
  } else stop("Choose survivorship or exponential exposures.")
  weight <- exp(log_weight - max(log_weight))
  data.frame(age = ages, baseline_exposure = weight / sum(weight), null_rate = h)
}

# --- "Assumptions and Fisher information" -----------------------------------
# Project the tested directions on the level/slope directions, using weighted
# least squares. QR avoids subtracting nearly equal information matrices.
# j is GM information for kappa; K adjusts BOTH kappa and gamma for level/slope.
# The frailty direction is evaluated at gamma = kappa = 0.
# Function: information_at_null
# Use the ages, exposure proportions and Gompertz parameters to calculate the
# information at kappa = gamma = 0. Return j for the Makeham-only test, K for
# the joint-boundary limit, their score correlation, and projection residuals.
# The zero-mass formula returned here applies to negative correlations only.
information_at_null <- function(age, baseline_exposure, lambda, beta) {
  h <- lambda * exp(beta * age)
  q <- lambda / beta * expm1(beta * age)
  root_weight <- sqrt(baseline_exposure / h)
  nuisance <- cbind(h, age * h) * root_weight
  tested <- cbind(kappa = rep(1, length(age)), gamma = -h * q) * root_weight
  residual <- qr.resid(qr(nuisance), tested)
  K <- crossprod(residual)
  rho <- K[1, 2] / sqrt(K[1, 1] * K[2, 2])
  stopifnot(all(eigen(K, symmetric = TRUE, only.values = TRUE)$values > 0))
  list(j = unname(K[1, 1]), K = K, rho = unname(rho),
       mass_zero = if (rho < 0) .5 + asin(rho) / (2 * pi) else NA_real_,
       weighted_residuals = residual)
}

# --- "Explicit form and critical values" -----------------------------------
# Columns of W are the adjusted scores for kappa and gamma, in that order.
# Under the null gamma is STILL estimated. Its edge maximum must be subtracted.
# We standardize K internally; positive rescaling preserves both constraints.
# Function: joint_statistic
# For each row of Gaussian scores W, compare the maximum over the nonnegative
# quadrant with the maximum on the gamma edge. Return one limiting statistic
# per draw. An interior candidate counts only when both coordinates are valid.
joint_statistic <- function(W, K) {
  scale <- sqrt(diag(K))
  C <- K / outer(scale, scale)
  w <- sweep(as.matrix(W), 2, scale, "/")
  z <- w %*% solve(C)
  edge_kappa <- pmax(w[, 1], 0)^2
  edge_gamma <- pmax(w[, 2], 0)^2
  interior <- rowSums(w * z)
  interior[rowSums(z >= 0) != 2] <- -Inf
  pmax(0, pmax(edge_kappa, edge_gamma, interior) - edge_gamma)
}

# --- Likelihood fitting: "Testing for a Makeham term" ------------------------
# The scientific model always uses the manuscript's age origin:
# mu(x) = kappa + lambda*exp(beta*x) /
#                   (1 + gamma*lambda/beta * (exp(beta*x)-1)).
# For numerical stability only, write t = x-min(x) and a=lambda*exp(beta*min(x)).
# We retain the original denominator; centering does NOT reset the frailty law.
# Internal parameters are log(a), log(beta*span), kappa/k_scale, gamma/g_scale.
# Zero kappa and gamma remain attainable; we do not replace them by tiny positives.
# Function: make_likelihood
# Prepare the rate, objective, gradient and convergence-check functions for
# one simulated death-count vector. All four fitted models use this same
# likelihood; their difference is which parameters are free to vary.
# Return these functions together with starting values and numerical bounds.
make_likelihood <- function(deaths, exposure, age, controls) {
  stopifnot(all(is.finite(deaths)), all(deaths >= 0),
            all(deaths == round(deaths)), all(is.finite(exposure)),
            all(exposure > 0), sum(deaths) > 0)
  x0 <- min(age)
  t <- age - x0
  span <- max(t)
  # A log-rate regression supplies starting values only, never the final fit.
  init <- lm.fit(cbind(1, t), log((deaths + .5) / exposure))$coefficients
  beta_start <- min(.3, max(.01, unname(init[2])))
  a_start <- sum(deaths) / sum(exposure * exp(beta_start * t))
  k_scale <- max(a_start, 1e-7)
  g_scale <- 1 / max(a_start / beta_start *
                      (exp(beta_start * max(t)) - exp(-beta_start * x0)), 1e-6)
  initial <- c(log(a_start), log(beta_start * span), 0, 0)
  lower <- c(-25, log(.0001 * span), 0, 0)
  upper <- c(3, log(.6 * span), 1e5, 1e5)

  # Convert the four scaled coordinates back to mortality parameters and rates.
  evaluate <- function(p) {
    a <- exp(p[1]); beta <- exp(p[2]) / span
    kappa <- p[3] * k_scale; gamma <- p[4] * g_scale
    et <- exp(beta * t); e0 <- exp(-beta * x0)
    q <- a / beta * (et - e0)
    s <- 1 + gamma * q
    g <- a * et / s
    # Derivative of q with respect to log(beta), holding a fixed.
    dq <- a / beta * ((beta * t - 1) * et + (1 + beta * x0) * e0)
    derivative <- cbind(g / s, g * (beta * t - gamma * dq / s),
                        rep(k_scale, length(t)), -g_scale * g * q / s)
    list(rate = kappa + g, derivative = derivative,
         parameters = c(lambda = a * exp(-beta * x0), beta = beta,
                        kappa = kappa, gamma = gamma))
  }
  # Half Poisson deviance is negative log likelihood plus a data-only constant.
  # It yields the SAME LRT while avoiding cancellation between large log likelihoods.
  objective <- function(p) {
    mu <- exposure * evaluate(p)$rate
    positive <- deaths > 0
    ratio <- (mu[positive] - deaths[positive]) / deaths[positive]
    term <- mu[positive] - deaths[positive] +
      deaths[positive] * (log(deaths[positive]) - log(mu[positive]))
    close <- abs(ratio) < .5
    term[close] <- deaths[positive][close] * (ratio[close] - log1p(ratio[close]))
    sum(mu[!positive]) + sum(term)
  }
  # Analytic derivatives let the optimizer avoid numerical differencing.
  gradient <- function(p) {
    v <- evaluate(p)
    colSums(v$derivative * as.numeric(exposure - deaths / v$rate))
  }
  # Scale scores by their information so tolerances are comparable across fits.
  check_fit <- function(p, active) {
    v <- evaluate(p)
    score <- gradient(p)[active]
    boundary <- active %in% c(3, 4) & p[active] <= controls$boundary_tolerance
    score[boundary] <- pmin(score[boundary], 0)
    information <- colSums(v$derivative[, active, drop = FALSE]^2 *
                            as.numeric(exposure / v$rate))
    max(abs(score) / sqrt(pmax(information, .Machine$double.eps)))
  }
  list(evaluate = evaluate, objective = objective, gradient = gradient,
       check_fit = check_fit, initial = initial, lower = lower, upper = upper)
}

# Function: fit_model
# Fit a model from several starting values. `active` identifies the free
# coordinates: 1 = level, 2 = slope, 3 = Makeham, 4 = frailty.
# `embedded` supplies nested fits as additional candidates. Return the selected
# parameters and fitted rates, plus the diagnostics for every attempted start.
fit_model <- function(likelihood, active, starts, controls, embedded = list()) {
  candidates <- list()
  attempts <- list()
  for (i in seq_along(starts)) {
    base <- starts[[i]]
    unpack <- function(z) { p <- base; p[active] <- z; p }
    ans <- tryCatch(nlminb(base[active],
      objective = function(z) likelihood$objective(unpack(z)),
      gradient = function(z) likelihood$gradient(unpack(z))[active],
      lower = likelihood$lower[active], upper = likelihood$upper[active],
      control = list(iter.max = controls$max_iterations, eval.max = 3000,
                     rel.tol = 1e-11, x.tol = 1e-9)), error = identity)
    if (inherits(ans, "error")) {
      attempts[[i]] <- data.frame(start = i, objective = NA_real_, convergence = 99L,
        score = NA_real_, iterations = NA_integer_, message = conditionMessage(ans))
      next
    }
    # PORT occasionally reports singular convergence at an otherwise good
    # solution. Refine such a result with a second bounded optimizer; retain
    # both messages rather than relabelling an unsuccessful fit as converged.
    if (ans$convergence != 0L) {
      first_message <- ans$message
      retry <- tryCatch(optim(ans$par,
        fn = function(z) likelihood$objective(unpack(z)),
        gr = function(z) likelihood$gradient(unpack(z))[active],
        method = "L-BFGS-B", lower = likelihood$lower[active],
        upper = likelihood$upper[active],
        control = list(maxit = controls$max_iterations, factr = 1e5, pgtol = 1e-8)),
        error = identity)
      if (!inherits(retry, "error") && retry$convergence == 0L &&
          retry$value <= ans$objective + controls$likelihood_tolerance) {
        ans$par <- retry$par
        ans$objective <- retry$value
        ans$convergence <- retry$convergence
        ans$message <- paste(first_message, "-> L-BFGS-B:", retry$message)
      }
    }
    p <- unpack(ans$par)
    score <- likelihood$check_fit(p, active)
    # Finite upper/lower safeguards are not scientific parameter restrictions.
    # If they determine a fit, flag it instead of reporting an ordinary result.
    safeguard <- any(p[active] >= likelihood$upper[active] - 1e-6) ||
      any(p[1:2] <= likelihood$lower[1:2] + 1e-6)
    ok <- ans$convergence == 0L && is.finite(ans$objective) &&
      is.finite(score) && score <= controls$score_tolerance && !safeguard
    attempts[[i]] <- data.frame(start = i, objective = ans$objective,
      convergence = ans$convergence, score = score, iterations = ans$iterations,
      message = paste(ans$message, if (safeguard) "[safeguard reached]" else ""))
    candidates[[length(candidates) + 1L]] <- list(p = p, objective = ans$objective,
      ok = ok, score = score, source = paste0("start_", i))
  }
  # Explicitly include nested fits, but only accept them as full optima if their
  # projected scores also satisfy the larger model's first-order conditions.
  for (fit in embedded) if (!is.null(fit$p)) {
    score <- likelihood$check_fit(fit$p, active)
    candidates[[length(candidates) + 1L]] <- list(p = fit$p,
      objective = fit$objective, ok = isTRUE(fit$ok) && score <= controls$score_tolerance,
      score = score, source = "nested_fit")
  }
  if (!length(candidates)) return(list(ok = FALSE, message = "All starts failed",
                                         attempts = do.call(rbind, attempts)))
  # Never silently select a poorer converged solution over a clearly better one.
  values <- vapply(candidates, function(z) z$objective, numeric(1))
  valid <- which(vapply(candidates, function(z) isTRUE(z$ok), logical(1)))
  best <- if (length(valid)) valid[which.min(values[valid])] else which.min(values)
  chosen <- candidates[[best]]
  if (chosen$objective > min(values) + controls$likelihood_tolerance) chosen$ok <- FALSE
  v <- likelihood$evaluate(chosen$p)
  chosen$parameters <- v$parameters
  chosen$rate <- v$rate
  chosen$attempts <- do.call(rbind, attempts)
  chosen$message <- if (chosen$ok) "ok" else "Check optimizer convergence, score or safeguard"
  chosen
}

# Function: likelihood_ratio
# Compare two nested fits using twice the difference in half-deviance.
# Return the raw statistic, the statistic after numerical-zero handling, and
# a status. Failed fits remain missing; they are not counted as non-rejections.
likelihood_ratio <- function(null, alternative, controls) {
  if (!isTRUE(null$ok) || !isTRUE(alternative$ok))
    return(list(T = NA_real_, raw_T = NA_real_, status = "fit_failed"))
  raw <- 2 * (null$objective - alternative$objective)
  if (raw < -controls$likelihood_tolerance)
    return(list(T = NA_real_, raw_T = raw, status = "nesting_failed"))
  # A finite optimizer cannot identify an exact atom. Record raw_T as well and
  # use a stated numerical tolerance for the empirical mass near zero.
  list(T = if (abs(raw) <= controls$zero_tolerance) 0 else max(0, raw),
       raw_T = raw, status = "ok")
}

# Gompertz alone has a particularly simple fit: for a fixed slope, the level
# has a closed-form estimate. The remaining score matches mean age at death
# to its fitted value and is monotone in beta. Solve it directly; this avoids
# iterative optimizer warnings at an already accurate two-parameter optimum.
# Function: fit_gompertz
# Fit the null Gompertz model by solving its slope score after profiling out
# the level. Return the same fit structure used by the other model fits.
fit_gompertz <- function(L, deaths, exposure, age, controls) {
  t <- age - min(age)
  target <- sum(t * deaths) / sum(deaths)
  slope_score <- function(beta) {
    w <- exposure * exp(beta * t)
    sum(t * w) / sum(w) - target
  }
  bounds <- exp(c(L$lower[2], L$upper[2])) / max(t)
  if (slope_score(bounds[1]) >= 0 || slope_score(bounds[2]) <= 0)
    stop("Gompertz slope is outside the numerical safeguards")
  beta <- uniroot(slope_score, bounds, tol = 1e-12)$root
  a <- sum(deaths) / sum(exposure * exp(beta * t))
  p <- c(log(a), log(beta * max(t)), 0, 0)
  v <- L$evaluate(p)
  score <- L$check_fit(p, 1:2)
  value <- L$objective(p)
  list(p = p, objective = value, ok = score <= controls$score_tolerance &&
         p[1] > L$lower[1] && p[1] < L$upper[1], score = score,
       source = "profile_level_and_slope_root", parameters = v$parameters,
       rate = v$rate, message = "Gompertz profile score solved",
       attempts = data.frame(start = 1L, objective = value, convergence = 0L,
         score = score, iterations = NA_integer_, message = "Monotone slope score root"))
}

# Function: fit_models
# Fit Gompertz (G) and Gompertz-Makeham (GM) to one simulated dataset.
# With joint = TRUE, also fit gamma-Gompertz (GG) and gamma-Gompertz-Makeham
# (GGM). The joint-boundary comparison is GG versus GGM: gamma is estimated
# on both sides, while the tested parameter is kappa alone.
fit_models <- function(deaths, exposure, age, joint = FALSE, controls) {
  L <- make_likelihood(deaths, exposure, age, controls)
  gompertz <- fit_gompertz(L, deaths, exposure, age, controls)
  p <- gompertz$p
  gm_start <- p; gm_start[3] <- .5; gm_start[2] <- gm_start[2] + log(1.15)
  gm <- fit_model(L, 1:3, list(p, gm_start), controls, list(gompertz))
  result <- list(G = gompertz, GM = gm,
                 GM_test = likelihood_ratio(gompertz, gm, controls))
  if (joint) {
    gg_start <- p; gg_start[4] <- .5; gg_start[2] <- gg_start[2] + log(1.15)
    gg <- fit_model(L, c(1, 2, 4), list(p, gg_start), controls, list(gompertz))
    starts <- list(p, gm$p, gg$p)
    extra <- p; extra[3:4] <- c(.5, .5); extra[2] <- extra[2] + log(1.25)
    starts <- Filter(Negate(is.null), c(starts, list(extra)))
    ggm <- fit_model(L, 1:4, starts, controls, list(gm, gg))
    result$GG <- gg; result$GGM <- ggm
    # IMPORTANT: GG is the null, not G. This tests kappa alone with gamma free.
    result$joint_test <- likelihood_ratio(gg, ggm, controls)
  }
  for (name in intersect(names(result), c("G", "GM", "GG", "GGM"))) {
    if (!is.null(result[[name]]$rate))
      result[[name]]$log_likelihood <- sum(dpois(deaths,
        exposure * result[[name]]$rate, log = TRUE))
  }
  result
}

# --- Saved summaries: keep failures visible ---------------------------------
# The conditional estimate excludes failed fits, so also save lower/upper
# bounds treating every failed fit as a non-rejection or a rejection.
# Function: probability_summary
# Summarize a logical event vector (for example, T > cutoff). NA records a
# failed test. Return its estimated probability, Monte Carlo standard error,
# Wilson interval, failure count, and bounds accounting for missing outcomes.
probability_summary <- function(event, total = length(event)) {
  good <- !is.na(event)
  n <- sum(good); yes <- sum(event[good])
  p <- if (n) yes / n else NA_real_
  # Wilson interval describes Monte Carlo uncertainty, not data uncertainty.
  z <- qnorm(.975)
  centre <- if (n) (p + z^2/(2*n))/(1 + z^2/n) else NA_real_
  half <- if (n) z*sqrt(p*(1-p)/n + z^2/(4*n^2))/(1 + z^2/n) else NA_real_
  data.frame(total = total, successful = n, failed = total - n, estimate = p,
    mcse = if (n) sqrt(p*(1-p)/n) else NA_real_,
    mc_lower = centre - half, mc_upper = centre + half,
    failure_lower = yes / total, failure_upper = (yes + total - n) / total)
}

# One checkpoint holds several independent replicates, including counts and fits.
# Its seed belongs to the checkpoint, not the worker: scheduling changes do not
# change results. B, chunk size and source-code versions are frozen on resume.
# Function: run_chunk
# Simulate and fit one batch for one (n, delta) scenario. Save the simulated
# deaths and all fit diagnostics in a checkpoint, then return a small progress
# summary. Writing a temporary file first prevents a partial save from looking
# like a completed checkpoint if R stops during the write.
run_chunk <- function(job, settings, design, scenarios, folder) {
  scenario <- scenarios[job$scenario, ]
  set.seed(job$seed)
  exposure <- scenario$n * design$baseline_exposure
  mean_deaths <- exposure * (design$null_rate + scenario$kappa)
  records <- vector("list", job$last - job$first + 1L)
  rows <- vector("list", length(records))
  for (i in seq_along(records)) {
    deaths <- rpois(nrow(design), mean_deaths)
    result <- tryCatch(fit_models(deaths, exposure, design$age,
      joint = scenario$delta == 0, controls = settings$controls), error = identity)
    failed <- inherits(result, "error")
    test <- if (failed) list(T = NA_real_, raw_T = NA_real_, status = conditionMessage(result)) else result$GM_test
    jt <- if (failed || scenario$delta != 0) list(T = NA_real_, raw_T = NA_real_, status = if (failed) "fit_failed" else "not_run") else result$joint_test
    rows[[i]] <- data.frame(scenario = job$scenario, replicate = job$first + i - 1L,
      n = scenario$n, delta = scenario$delta, kappa = scenario$kappa,
      is_threshold = scenario$is_threshold, total_deaths = sum(deaths),
      GM_T = test$T, GM_raw_T = test$raw_T, GM_status = test$status,
      joint_T = jt$T, joint_raw_T = jt$raw_T, joint_status = jt$status)
    records[[i]] <- list(deaths = deaths, fits = if (failed) NULL else result,
                        error = if (failed) conditionMessage(result) else NULL)
  }
  out <- list(job = job, scenario = scenario, exposure = exposure,
              generating_means = mean_deaths, results = do.call(rbind, rows),
              records = records)
  file <- file.path(folder, sprintf("chunk_%05d.rds", job$id))
  temporary <- paste0(file, ".tmp")
  saveRDS(out, temporary, compress = TRUE)
  if (!file.rename(temporary, file)) stop("Could not finalize checkpoint: ", file)
  data.frame(job = job$id, replicates = length(records),
             GM_failed = sum(out$results$GM_status != "ok"),
             joint_failed = sum(!out$results$joint_status %in% c("ok", "not_run")))
}

# Function: results_folder
# Use results/ directly. dir.create also works when the folder already exists;
# showWarnings = FALSE suppresses that harmless message.
# Keep the lock check so table/figure scripts do not read an unfinished run.
results_folder <- function(project_dir) {
  folder <- file.path(project_dir, "results")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(file.path(folder, ".running")))
    stop("The results folder is locked by a simulation. Wait for it to finish before reading or changing outputs.")
  folder
}
