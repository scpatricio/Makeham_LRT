# When does a mortality schedule need a Makeham term?

Code accompanying the research note *When does a mortality schedule need a Makeham term? The standard test answers at half the stated level* (Silvio C. Patricio, Interdisciplinary Center on Population Dynamics, University of Southern Denmark).

## About the note

The Makeham constant cannot be negative, so testing whether it is zero places the parameter on the boundary of its range. Under independent Poisson death counts, the likelihood-ratio statistic then converges to an equal mixture of a point mass at zero and a chi-squared distribution with one degree of freedom, rather than to the usual chi-squared limit. A test referred to the chi-squared distribution with one degree of freedom therefore runs at half its nominal level. The note derives this result for a general additive hazard and verifies its conditions for the Gompertz–Makeham and gamma-Gompertz–Makeham models.

This repository contains the code for the Monte Carlo study reported in the note, which checks the limiting distribution and the size loss from the conventional calibration.

## Simulation design

Data sets are simulated under the null hypothesis from a Gompertz model with no Makeham term:

- single-year ages 30 to 89
- λ = 10⁻⁵, β = 0.10, κ = 0
- baseline exposures proportional to the implied survivorship, normalized to one at age 30 and multiplied by 10⁵ (about 5.4 million person-years and 53,600 expected deaths in total)
- death counts drawn as independent Poisson variables
- 20,000 replicates

For each replicate, the Gompertz model is fitted by profiling λ out in closed form and maximizing over β, the Gompertz–Makeham model is fitted over κ ≥ 0 from several starting values, and the likelihood-ratio statistic is computed. Values below 10⁻⁶ are treated as zero.

The restricted fit needs to be accurate. If the Gompertz maximum is found only approximately, cases that should give a statistic of exactly zero return small positive values instead, and the point mass at zero is underestimated.

## Results

| | Value |
|---|---|
| Replicates with statistic equal to zero | 50.2% |
| Rejection rate, conventional critical value 3.8415 | 2.5% |
| Rejection rate, corrected critical value 2.7055 | 5.2% |

Conditional on being positive, the quantiles of the statistic track those of the chi-squared distribution with one degree of freedom.

## Reproducing the results

[Add the software version, required packages, and the command that runs the simulation.]

## Citation

[Add the citation once the note is published.]
