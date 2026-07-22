# cran-comments.md — xtife 0.1.5

## Submission note

This is a fast follow-up to 0.1.4 (submitted earlier today). It aligns two
implementation details of the unbalanced-panel estimator with the reference
methodology of Su, Wang and Wang (2025) <doi:10.2139/ssrn.5177283> and fixes
a related defect:

* The nuclear-norm penalty grid now uses c * max(N, T) (the grid of Su,
  Wang and Wang 2025), replacing c * sqrt(max(N, T)).
* The default initialisation of ife_unbalanced() is now the
  nuclear-norm-regularised initial estimator (init = "nnr"), the first step
  of the reference two-step procedure.
* Fixed a defect where a penalty exceeding the largest singular value
  shrank the low-rank component to exactly zero, making the singular-value-
  thresholding factor-selection rule degenerate (it could return
  r_hat = min(N, T)). Degenerate penalty candidates are now excluded.

Apologies for the quick resubmission; the fix affects default behaviour,
so we preferred to correct it before the package enters wide use.

## Test environments

* local: macOS (Apple Silicon, aarch64-apple-darwin20), R 4.5.2

## R CMD check results

0 errors | 0 warnings | 1 note

* NOTE: "Days since last update: 0" — see submission note above.

Any additional URL/DOI NOTEs observed locally were transient network
timeouts (libcurl error 28) on stable URLs (https://www.gnu.org/licenses/
and https://doi.org/), not genuinely invalid links.

## Downstream dependencies

None.
