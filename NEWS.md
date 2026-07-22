# xtife 0.1.5

Alignment of the unbalanced-panel estimator with the reference methodology
of Su, Wang and Wang (2025) <doi:10.2139/ssrn.5177283>, plus a related fix.

* The nuclear-norm penalty grid is now `c * max(N, T)` for
  `c in {0.01, 0.1, 1, 10}` — the grid of Su, Wang and Wang (2025,
  footnote 12). Previously `c * sqrt(max(N, T))`.
* The default initialisation of `ife_unbalanced()` is now `init = "nnr"`
  (the nuclear-norm-regularised consistent initial estimator, the first
  step of the reference two-step procedure). `init = "ols"` remains
  available as a faster warm start.
* Fixed: a penalty exceeding the largest singular value shrank the
  low-rank component to exactly zero, making the singular-value-
  thresholding rule in `ife_select_r_unb()` degenerate (it could return
  `r_hat = min(N, T)`). Degenerate penalty candidates are now excluded,
  an all-degenerate grid errors informatively, and the SVT count guards
  against a zero first singular value.

# xtife 0.1.4

* `ife_unbalanced()` gains full additive fixed-effect support via `force`
  (`"none"`, `"unit"`, `"time"`, `"two-way"`), estimated jointly with the
  factors by EM on the imputed panel; replaces the previous grand-mean
  centering, which was selection-biased under informative missingness.
* Fixed the sign convention of the static bias correction (Bai 2009:
  corrected = raw − B/N − C/T) and the premultiplication of the
  Moon–Weidner dynamic bias-correction contributions.
* Standard-error degrees of freedom account for additive fixed effects
  under `force`; factor space augmented with the additive effects for
  SE and bias-correction computations.
* Four standard-error types for unbalanced panels: `"standard"`,
  `"robust"` (HC1), `"cluster"` (by unit), `"hac"` (Bartlett kernel,
  bandwidth `floor(2 * T^0.2)`).
* Documentation cites Su, Wang and Wang (2025) in standard author–year
  form throughout.
