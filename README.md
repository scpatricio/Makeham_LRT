# Does a mortality schedule need a Makeham term? Calibrating the likelihood-ratio test

R code for the numerical illustrations in **Does a mortality schedule need a Makeham term? Calibrating the likelihood-ratio test**, by Silvio C. Patricio.

The simulations illustrate three questions: how the Makeham likelihood-ratio test behaves when the constant is absent, how often it detects a positive constant, and how its null distribution changes when a frailty variance is estimated at a true value of zero. All death counts are simulated; no external mortality data are needed.

## Files and requirements

| File | What it does |
|---|---|
| `01_functions.R` | Defines mortality rates, likelihood fitting, information calculations, the joint-boundary statistic and checkpoint functions. |
| `02_run_simulations.R` | Sets the design, runs simulations in parallel and saves counts, fits and diagnostics. |
| `03_make_tables.R` | Reads saved results and creates the calibration table and supporting numerical summaries. |
| `04_make_figures.R` | Reads saved results and creates the two manuscript figures with ggplot2. |

The calculations use base R and the bundled `parallel` package. Only the figures require an additional package:

```r
install.packages("ggplot2")  # Run once, if not already installed.
```

The figure script uses `linewidth`, available in ggplot2 3.4.0 and later. The scripts record the R environment in `results/session.txt`; save that file with results used for publication.

## Run the illustrations

Set R's working directory to this repository folder. In RStudio, creating a project in this folder is a convenient way to do that. Then run:

```r
Sys.setenv(MAKEHAM_MODE = "full", MAKEHAM_WORKERS = "8")
source("02_run_simulations.R")
source("03_make_tables.R")
source("04_make_figures.R")
```

`02_run_simulations.R` sources the functions automatically. Tables and figures use the saved simulation settings, so they can be regenerated in a new R session without fitting models again.

For a short trial, set `MAKEHAM_MODE="smoke"` before running the same sequence. The modes use:

| Mode | Datasets per scenario | Gaussian draws for the joint limit |
|---|---:|---:|
| `smoke` | 20 | 5,000 |
| `pilot` | 200 | 100,000 |
| `full` (default) | 10,000 | 100,000 |

Use a separate results folder for each design or mode: after a run finishes, move or rename `results/` before starting a different run. Tables and figures always read the active `results/` folder.

## Manuscript reference

The reference for this repository is `Makeham_note-7.pdf`. Section 3.1 (page 7) specifies ages **30 through 90 inclusive** and **10,000 independent simulated datasets per exposure level** for calibration. Section 3.3 (page 9) specifies **100,000 Gaussian draws** for the joint-boundary limit. These Gaussian draws are not additional simulated mortality datasets.

The manuscript settings correspond to `MAKEHAM_MODE="full"`. Smoke and pilot modes are optional short checks and do not reproduce the manuscript's replication count.

| Manuscript item | Script location |
|---|---|
| Section 3.1: Gompertz null, ages, exposures and 10,000 samples | `02_run_simulations.R`: `settings`; `01_functions.R`: `make_design`, `run_chunk` |
| Section 3.2: local alternatives and 80% power threshold | `02_run_simulations.R`: `threshold_delta`, `scenarios$kappa` |
| Section 3.3: estimated frailty under both hypotheses | `01_functions.R`: `fit_models`, GG versus GGM |
| Section 3.3: 100,000 Gaussian draws using the true information matrix | `02_run_simulations.R`: `B_limit`, `W`, `T0` |
| Table 1 | `03_make_tables.R`: `table1_calibration.tex` |
| Figure 1: empirical and local asymptotic power | `04_make_figures.R`: `figure_power` |
| Figure 2: joint-boundary upper-tail probabilities | `04_make_figures.R`: `figure_joint` |

## Simulation design

The settings near the top of `02_run_simulations.R` reproduce the numerical illustration design:

- Single-year ages 30–90, with Gompertz parameters lambda = 0.00002 and beta = 0.10.
- Exposure proportions based on survivorship under that Gompertz schedule, normalized to sum to one.
- Total exposure n = 10^4, 10^5, 10^6 and 10^7 person-years.
- Standardized shifts delta = 0, 0.25, ..., 3, plus the local 80% power threshold (approximately 2.486).
- Test level alpha = 0.05 and target power = 0.80.
- Master seed 1; batches of 100 datasets in full/pilot mode and 10 in smoke mode.

The supplied code uses 10,000 datasets at each (n, delta) combination in full mode. There are 56 scenarios: four exposure totals and 14 shifts. The full run gives 560,000 Gompertz versus Gompertz–Makeham comparisons. The 40,000 samples at delta = 0 also supply the calibration and joint-boundary comparisons; these reuse the same counts rather than generating new null datasets.

The positive Makeham constant is `kappa = delta / sqrt(n * j_eff)`. Thus, at a fixed standardized shift, the constant decreases as exposure grows. The information calculation gives approximately `j_eff = 243.9637` and adjusted-score correlation `rho = -0.6301555` for this design, matching the rounded values in Sections 3.2 and 3.3. The grid spacing, seed and batch size are implementation settings from the supplied scripts; the manuscript does not specify them.

Comments in the scripts refer to manuscript section titles rather than equation numbers.

## Parallel execution and resuming

The simulations run in independent PSOCK R sessions, including when called from RStudio on macOS. By default the supplied runner uses the detected physical core count minus two, with a minimum of one. Set `MAKEHAM_WORKERS` explicitly to control this; eight workers is a reasonable starting configuration on an Apple M4 Max with 36 GB of memory. It is not a measured optimum.

Each worker handles one batch at a time, saves its detailed checkpoint, and returns a small progress summary. The seed belongs to the batch, so changing the worker count does not change the simulated counts. Changing batch size does change seed assignment. Math-library threads are limited to one per worker to avoid competing layers of parallel work.

Rerun the simulation script with unchanged settings and source files to resume an interrupted run. Completed checkpoints are reused. The manifest compares settings, source-file hashes and R version before resuming. Comment edits also change a source hash; use a separate run after editing either simulation source file.

A `results/.running/` lock prevents two sessions from writing to the same run and stops post-processing while the run is active. If R crashes, remove this directory only after confirming the original simulation has stopped.

This repository uses `results/` directly. It does not migrate older nested folders or create a `code_snapshot` folder. Existing flat results can be read by the table and figure scripts without rerunning simulations. The rewritten source files will not resume a run created with the older source hashes; keep the original scripts with that run if resumption is needed.

## Outputs

```text
results/
  manifest.rds
  design.csv
  scenarios.csv
  joint_information_K.csv
  joint_limit.rds
  replicate_results.rds
  replicate_results.csv
  session.txt
  worker_progress.log          # Parallel runs
  checkpoint_timings.csv
  checkpoints/
    chunk_00001.rds
    ...
  tables/
    table1_calibration.tex
    table1_calibration_details.csv
    all_probability_summaries.csv
    all_probability_summaries.rds
    power_curves.csv
    detection_threshold.csv
    joint_boundary_summary.csv
    joint_limit_summary.csv
    numbers_for_manuscript.csv
    fitted_parameters.csv
    zero_tolerance_sensitivity.csv
  figures/
    figure1_local_power.pdf
    figure2_joint_boundary.pdf
```

The `N` in the figure legends denotes the total exposure called `n` in the text, as in the attached manuscript figures. The `.tex` table is a snippet to include in the manuscript. `numbers_for_manuscript.csv` contains information values and critical values; `detection_threshold.csv` contains power at the target threshold. The figure script preserves the supplied colours, shapes, legends and layout, and saves PDF files only. Monte Carlo error bars are available as a commented plotting layer and are off by default.

For a particular simulated dataset:

```r
chunk <- readRDS("results/checkpoints/chunk_00001.rds")
chunk$results[1, ]                    # Test statistics and statuses
chunk$records[[1]]$deaths             # Simulated age-specific death counts
chunk$exposure                       # Person-years at each age
chunk$generating_means                # Expected deaths for this scenario
chunk$records[[1]]$fits$GM$parameters  # Selected parameter estimates
chunk$records[[1]]$fits$GM$attempts    # Diagnostics for each starting point
```

GG and GGM fits are available only for the zero-shift scenarios. Checkpoints retain fitted rates, log likelihoods, optimizer messages and raw likelihood-ratio statistics as well as the selected estimates. Full-run files can be large. The included `.gitignore` excludes generated results and local R session files; the scripts recreate the output folders. Keep the original publication results separately if you want an exact record of that run.

## Reading the results

The optimizer minimizes half the Poisson deviance, which gives the same likelihood-ratio statistic as maximizing the Poisson log likelihood. Age centering improves numerical scaling while retaining the original gamma-frailty denominator and parameter interpretation. Both zero Makeham and zero frailty remain attainable.

The Makeham-only comparison fits G versus GM. The joint-boundary comparison fits **GG versus GGM**, estimating frailty under both hypotheses even though its true value is zero. The latter tests only the Makeham term; comparing G versus GGM would be a different test.

The joint critical value uses the **known generating information matrix K**. The experiment isolates the reference-distribution effect; it does not validate a rule that selects a reference from the fitted frailty estimate or establish finite-sample performance for an estimated-K cutoff.

Convergence, projected-score and nesting checks remain in the code because failed fits can distort estimated rejection probabilities. Multiple starting values help, but do not prove that every fit reaches the global optimum. Failed comparisons receive `NA`. Summaries report probabilities among successful fits, failure counts and bounds that treat failed comparisons as either rejections or non-rejections.

Monte Carlo standard errors and Wilson intervals describe simulation uncertainty. The local detection threshold is an asymptotic approximation, so finite-exposure power at that threshold need not equal 80%. In the joint-limit summary, the rejection probability at the simulated quantile uses the same Gaussian draws that defined it; that row is not an independent calibration check.