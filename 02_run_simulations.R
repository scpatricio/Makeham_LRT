setwd("~/Desktop/fertility/Makeham_test/numerical_illustrations/Makeham_LRT/")
rm(list = ls())

# Run this file first. The other scripts only read its saved results.
# In RStudio: open this folder as a project, then source this script.
# In Terminal: Rscript 02_run_simulations.R
# A quick check: MAKEHAM_MODE=smoke Rscript 02_run_simulations.R
# All paths are relative to this script when run through Rscript; when sourced,
# set the working directory to this folder (or use source(..., chdir = TRUE)).

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project_dir <- if (length(script_argument)) dirname(normalizePath(sub("^--file=", "", script_argument[1]))) else getwd()
source(file.path(project_dir, "01_functions.R"))

# --- Choices to edit ---------------------------------------------------------
mode <- Sys.getenv("MAKEHAM_MODE", "full")  # "smoke", "pilot", or "full"
if (!mode %in% c("smoke", "pilot", "full")) stop("Unknown mode")
B <- switch(mode, smoke = 20L, pilot = 200L, full = 10000L)
B_limit <- switch(mode, smoke = 5000L, pilot = 100000L, full = 100000L)

# --- Simulation design: "Finite-sample calibration" -------------------------
# Ages, hazards and exposure proportions follow the manuscript illustration.
# Change the values here for a separate experiment. The script recalculates
# information, local alternatives and critical values from these settings.
settings <- list(
  ages = 30:90,
  lambda = 2e-5,
  beta = .10,
  exposure_pattern = "survivorship",  # alternatively "exponential"
  exposure_decay = NULL,              # required for exponential exposures
  exposure_totals = c(1e4, 1e5, 1e6, 1e7),
  delta_grid = seq(0, 3, by = .25),
  alpha = .05,
  target_power = .80,
  B = B,
  B_limit = B_limit,
  seed = 1L,
  chunk_size = if (mode == "smoke") 10L else 100L,
  controls = list(max_iterations = 1000L, score_tolerance = 1e-4,
                  boundary_tolerance = 1e-9, likelihood_tolerance = 1e-6,
                  zero_tolerance = 1e-7)
)
# The full run fits hundreds of thousands of datasets. Start with pilot mode
# to assess runtime and failure rates on your own computer.
cores <- parallel::detectCores(logical = FALSE)
if (is.na(cores)) cores <- 2L
# Use physical cores minus two by default, preserving the supplied runner's
# choice. On the M4 Max, MAKEHAM_WORKERS=8 is a useful starting setting if you
# also want to use RStudio while simulations run. This changes scheduling only.
workers <- max(1L, cores - 2L)
workers <- as.numeric(Sys.getenv("MAKEHAM_WORKERS", workers))
stopifnot(is.finite(workers), workers >= 1L, workers == floor(workers))
workers <- as.integer(workers)
# One results directory per project. Keep a separate project copy if you want
# several designs or pilot/full runs; incompatible results are never overwritten.
output_dir <- results_folder(project_dir)

# Each worker receives the fixed design once, then only a small job descriptor.
# Results and simulated counts go straight to separate checkpoint files; only
# a small diagnostic row travels back to the main R process.
# Function: run_worker_chunk
# Give one checkpoint job to a worker, record its elapsed time, and append a
# short completion message to the worker log. Detailed fits stay on disk.
run_worker_chunk <- function(job) {
  context <- get("simulation_context", envir = .GlobalEnv)
  started <- proc.time()[3]
  report <- run_chunk(job, context$settings, context$design,
                      context$scenarios, context$folder)
  report$elapsed_seconds <- unname(proc.time()[3] - started)
  cat(format(Sys.time(), "%H:%M:%S"), " finished checkpoint ", job$id,
      "; GM failures=", report$GM_failed,
      "; joint failures=", report$joint_failed, "\n", sep = "")
  report
}

# Function: start_workers
# Open independent R sessions with PSOCK, which works on macOS and in RStudio.
# Limit each session to one math-library thread to avoid nested parallel work.
start_workers <- function(number, log_file) {
  # Avoid multiple processes each asking a threaded math library for many threads.
  # These settings are inherited at worker startup, then restored in this R
  # session. Matrix operations here are small, so one math thread is sufficient.
  variables <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS")
  old <- Sys.getenv(variables, unset = NA_character_)
  on.exit({
    for (name in variables) {
      if (is.na(old[[name]])) Sys.unsetenv(name)
      else do.call(Sys.setenv, setNames(list(old[[name]]), name))
    }
  }, add = TRUE)
  do.call(Sys.setenv, setNames(as.list(rep("1", length(variables))), variables))
  # PSOCK uses separate R sessions and works with RStudio on macOS.
  parallel::makePSOCKcluster(number, outfile = log_file)
}

# Function: run_simulations
# Build the design and scenario grid, save a manifest, run unfinished batches,
# and combine their compact result tables. The manifest and per-batch seeds
# allow an interrupted run to resume without changing its simulated datasets.
run_simulations <- function() {
  chunk_dir <- file.path(output_dir, "checkpoints")
  dir.create(chunk_dir, showWarnings = FALSE)
  # Prevent two R sessions writing into the same run simultaneously.
  lock <- file.path(output_dir, ".running")
  if (!dir.create(lock, showWarnings = FALSE))
    stop("This run is locked. If an earlier R session crashed and is no longer running, remove the .running folder.")
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)

  design <- make_design(settings$ages, settings$lambda, settings$beta,
                        settings$exposure_pattern, settings$exposure_decay)
  info <- information_at_null(design$age, design$baseline_exposure,
                               settings$lambda, settings$beta)
  threshold_delta <- qnorm(1 - settings$alpha) + qnorm(settings$target_power)
  delta <- sort(unique(c(settings$delta_grid, threshold_delta)))
  scenarios <- expand.grid(n = settings$exposure_totals, delta = delta)
  scenarios$scenario <- seq_len(nrow(scenarios))
  scenarios$kappa <- scenarios$delta / sqrt(scenarios$n * info$j)
  scenarios$is_threshold <- abs(scenarios$delta - threshold_delta) < 1e-12
  # delta=0 is shared by calibration, power and the joint-boundary experiment.
  # Sharing these null datasets saves fits and permits paired comparisons.
  scenarios$theory_boundary <- pnorm(scenarios$delta - qnorm(1 - settings$alpha))
  scenarios$theory_usual <- pnorm(scenarios$delta - qnorm(1 - settings$alpha/2))

  # Each job fits one batch at one exposure total and one standardized shift.
  # Job order also fixes seed assignment; keep it unchanged when reproducing a run.
  jobs <- list()
  for (s in seq_len(nrow(scenarios))) {
    for (first in seq.int(1L, settings$B, by = settings$chunk_size)) {
      jobs[[length(jobs) + 1L]] <- list(id = length(jobs) + 1L, scenario = s,
        first = first, last = min(settings$B, first + settings$chunk_size - 1L))
    }
  }
  # Assign seeds to jobs rather than workers, so worker count does not affect
  # the random counts. Reserve one additional seed for the Gaussian limit.
  set.seed(settings$seed)
  seeds <- sample.int(.Machine$integer.max, length(jobs) + 1L, replace = FALSE)
  for (j in seq_along(jobs)) jobs[[j]]$seed <- seeds[j]
  # Retain a simple resume check. Even a comment edit changes the file hash;
  # it is safer to start a separate run than combine files from different versions.
  code <- file.path(project_dir, c("01_functions.R", "02_run_simulations.R"))
  signature <- list(settings = settings, code_md5 = unname(tools::md5sum(code)),
                    R_version = R.version.string)
  manifest_file <- file.path(output_dir, "manifest.rds")
  if (file.exists(manifest_file)) {
    previous <- readRDS(manifest_file)
    if (!identical(signature, previous$signature)) {
      stop("Settings, source files or R version changed. Move the existing results folder aside or use a separate project folder; old results were not overwritten.")
    }
  } else {
    saveRDS(list(signature = signature, design = design, information = info,
                 scenarios = scenarios, jobs = jobs, limit_seed = tail(seeds, 1),
                 created = Sys.time()), manifest_file)

  }
  write.csv(design, file.path(output_dir, "design.csv"), row.names = FALSE)
  write.csv(scenarios, file.path(output_dir, "scenarios.csv"), row.names = FALSE)
  write.csv(info$K, file.path(output_dir, "joint_information_K.csv"))
  writeLines(c(capture.output(sessionInfo()), paste("Workers:", workers),
               paste("Mode:", mode)), file.path(output_dir, "session.txt"))

  # --- "The joint boundary at zero frailty": simulate its Gaussian limit ---
  # c_K uses the TRUE design information. This isolates finite-exposure
  # approximation error. It is not a study of an estimated-K calibration rule.
  limit_file <- file.path(output_dir, "joint_limit.rds")
  if (!file.exists(limit_file)) {
    set.seed(tail(seeds, 1))
    W <- matrix(rnorm(2 * settings$B_limit), ncol = 2) %*% chol(info$K)
    T0 <- joint_statistic(W, info$K)
    saveRDS(list(W = W, T0 = T0, K = info$K,
                critical = unname(quantile(T0, 1 - settings$alpha, type = 1)),
                seed = tail(seeds, 1)), limit_file)
  }
  limit <- readRDS(limit_file)
  message("Information j = ", signif(info$j, 7), "; rho = ", signif(info$rho, 5),
          "; joint critical value = ", signif(limit$critical, 5))

  # Check completed files before resuming. A damaged checkpoint stops the run
  # explicitly; it must not be silently mistaken for completed simulations.
  done <- vapply(jobs, function(job) {
    path <- file.path(chunk_dir, sprintf("chunk_%05d.rds", job$id))
    if (!file.exists(path)) return(FALSE)
    saved <- readRDS(path)
    if (!identical(saved$job, job)) stop("Checkpoint does not match job: ", path)
    TRUE
  }, logical(1))
  todo <- jobs[!done]
  message(sum(done), "/", length(jobs), " checkpoints already complete.")
  if (length(todo)) {
    active_workers <- min(workers, length(todo))
    started <- proc.time()[3]
    message("Running ", length(todo), " checkpoints using ", active_workers,
            " worker(s). Each checkpoint has up to ", settings$chunk_size,
            " independent replicates.")
    if (active_workers == 1L) {
      report <- lapply(todo, function(job) {
        begin <- proc.time()[3]
        row <- run_chunk(job, settings, design, scenarios, chunk_dir)
        row$elapsed_seconds <- unname(proc.time()[3] - begin)
        message("Finished checkpoint ", job$id, "/", length(jobs))
        row
      })
    } else {
      log_file <- file.path(output_dir, "worker_progress.log")
      message("Worker progress is written to: ", log_file)
      cluster <- start_workers(active_workers, log_file)
      on.exit(parallel::stopCluster(cluster), add = TRUE)
      parallel::clusterCall(cluster, function(path, context) {
        source(path, local = .GlobalEnv)
        assign("simulation_context", context, envir = .GlobalEnv)
        NULL
      }, file.path(project_dir, "01_functions.R"),
      list(settings = settings, design = design, scenarios = scenarios, folder = chunk_dir))
      # One continuously load-balanced queue: no barriers between waves of jobs.
      # chunk.size=1 means one checkpoint per assignment, NOT one replicate.
      # run_chunk sets its own seed, so worker count and ordering do not affect
      # the generated samples. No nested parallelism inside an optimizer.
      report <- parallel::parLapplyLB(cluster, todo, run_worker_chunk, chunk.size = 1L)
    }
    report <- do.call(rbind, report)
    seconds <- unname(proc.time()[3] - started)
    # Append session timing so a resumed run retains the earlier timing rows.
    report$workers <- active_workers
    report$session_started <- format(Sys.time() - seconds, "%Y-%m-%d %H:%M:%S")
    timing_file <- file.path(output_dir, "checkpoint_timings.csv")
    already_exists <- file.exists(timing_file)
    write.table(report, timing_file, sep = ",", row.names = FALSE,
                col.names = !already_exists, append = already_exists)
    message("Finished ", length(todo), " checkpoints in ", round(seconds/60, 2),
            " minutes; failed tests: GM=", sum(report$GM_failed),
            ", joint=", sum(report$joint_failed))
  }
  # Read only one large checkpoint at a time; retain its compact results table.
  rows <- lapply(jobs, function(job) {
    readRDS(file.path(chunk_dir, sprintf("chunk_%05d.rds", job$id)))$results
  })
  results <- do.call(rbind, rows)
  results$GM_p_half <- ifelse(is.na(results$GM_T), NA_real_,
    ifelse(results$GM_T == 0, 1, .5 * pchisq(results$GM_T, 1, lower.tail = FALSE)))
  results$GM_reject_boundary <- results$GM_T > qchisq(1 - 2 * settings$alpha, 1)
  results$GM_reject_usual <- results$GM_T > qchisq(1 - settings$alpha, 1)
  results$joint_reject_boundary <- results$joint_T > qchisq(1 - 2 * settings$alpha, 1)
  results$joint_reject_design_critical <- results$joint_T > limit$critical
  saveRDS(results, file.path(output_dir, "replicate_results.rds"))
  write.csv(results, file.path(output_dir, "replicate_results.csv"), row.names = FALSE)
  message("Saved ", nrow(results), " GM tests and ", sum(results$delta == 0),
          " joint-boundary tests to: ", output_dir)
  message("Next: run 03_make_tables.R, then 04_make_figures.R. Both read directly from results/.")
  invisible(results)
}

simulation_results <- run_simulations()
