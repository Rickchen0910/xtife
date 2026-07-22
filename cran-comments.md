# cran-comments.md — xtife 0.1.4

## Test environments

* local: macOS (Apple Silicon, aarch64-apple-darwin20), R 4.5.2

## R CMD check results

0 errors | 0 warnings | 1 note

* NOTE: "New submission" (this is the first submission of xtife to CRAN).

Any additional URL/DOI NOTEs observed locally were transient network
timeouts (libcurl error 28) on stable URLs (https://www.gnu.org/licenses/
and https://doi.org/), not genuinely invalid links.

## Package overview

xtife implements the interactive fixed effects panel estimator of
Bai (2009) <doi:10.3982/ECTA6135> for balanced and unbalanced panels,
with analytical standard errors, asymptotic bias correction, and factor
number selection. Unbalanced-panel estimation and inference follow
Su, Wang and Wang (2025) <doi:10.2139/ssrn.5177283>. The package uses
base R only (Imports: stats) with no compiled code.

## Downstream dependencies

None — this is a new package.
