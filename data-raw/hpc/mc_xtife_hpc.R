# ==============================================================================
# mc_xtife_hpc.R
# HPC worker script: runs ONE (dgp_id, N, TT, fill) parameter cell
# and saves results to data-raw/hpc/results/.
#
# Usage (called by run_mc_job.sh):
#   Rscript mc_xtife_hpc.R <dgp_id> <N> <TT> <fillx100> <ncores> [B]
#
#   dgp_id  : 1 = static balanced,  2 = static unbalanced,
#             3 = dynamic balanced,  4 = dynamic unbalanced
#   N       : number of cross-sectional units
#   TT      : number of time periods
#   fillx100: e.g. 80 means 80 % of cells observed (used for DGPs 2 & 4)
#   ncores  : parallel workers ($NSLOTS from SGE)
#   B       : replications (default 1000)
#
# Output: $XTIFE_MC_BASE/results/dgp{d}_N{N}_T{TT}_f{fill}_{B}reps.rds
#   List: params, summ (summary stats), raw_mat (Bx21 matrix), timing
# ==============================================================================

# Load xtife -------------------------------------------------------------------
# Strategy: try the installed package first; if unavailable (no system-wide
# install on HPC), source the two R files directly from xtife_src/.
# To use the source approach, upload R/ife.R and R/ife_unbalanced.R to:
#   /home/bc25911/xtife_mc/xtife_src/
# No installation, no admin rights, no internet access needed on compute nodes.
# ------------------------------------------------------------------------------
.xtife_base <- Sys.getenv("XTIFE_MC_BASE", unset = "/home/bc25911/xtife_mc")
.xtife_src  <- .xtife_base   # ife.R and ife_unbalanced.R sit directly in xtife_mc/

if (requireNamespace("xtife", quietly = TRUE)) {
  suppressPackageStartupMessages(library(xtife))
  cat("[xtife] loaded from installed package\n")
} else if (dir.exists(.xtife_src)) {
  for (.f in file.path(.xtife_src, c("ife.R", "ife_unbalanced.R")))
    source(.f, local = FALSE)
  cat(sprintf("[xtife] loaded from source: %s\n", .xtife_src))
} else {
  stop(
    "Cannot load xtife.\n",
    "Option A - install to user library on HPC:\n",
    "  mkdir -p ~/R/library\n",
    "  R -e 'install.packages(\"xtife\", lib=\"~/R/library\")'\n",
    "  Then add: export R_LIBS_USER=~/R/library to run_mc_job.sh\n\n",
    "Option B - source from files (no installation needed):\n",
    "  rsync -av R/ife.R R/ife_unbalanced.R ",
    "bc25911@hpc.essex.ac.uk:", .xtife_src, "/\n"
  )
}
rm(.xtife_base, .xtife_src, list = ls(pattern = "^\\.f$"))

# ==============================================================================
# 0 - Command-line arguments & fixed parameters
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L)
  stop("Usage: Rscript mc_xtife_hpc.R <dgp_id> <N> <TT> <fillx100> <ncores> [B]")

dgp_id   <- as.integer(args[1L])
N        <- as.integer(args[2L])
TT       <- as.integer(args[3L])
fill_int <- as.integer(args[4L])    # e.g. 80
nc       <- as.integer(args[5L])    # $NSLOTS
B        <- if (length(args) >= 6L) as.integer(args[6L]) else 1000L
fill     <- fill_int / 100          # e.g. 0.80

## Fixed across all cells
r_true <- 2L
r_max  <- 5L

## DGP-specific beta_true:
##   Static  DGPs (1, 2): beta = 1.0.
##   Dynamic DGPs (3, 4): beta = 0.5 (|beta| < 1 for stationarity, required by
##                        Moon & Weidner 2017; beta = 1 would be a unit root).
beta_true <- if (dgp_id %in% c(1L, 2L)) 1.0 else 0.5

## Error structure (env var MC_ERR), four variants:
##   "homo"      (default) i.i.d. homoskedastic, serially uncorrelated.
##   "hetero"    factor-driven heteroskedasticity (sd tied to |lambda_i'F_t|),
##               serially uncorrelated.
##   "serial"    homoskedastic but serially correlated MA(1) errors.
##   "serialhet" heteroskedastic AND serially correlated (SWW2025 DGP 1).
## Heteroskedasticity correlates the error variance with the regressor ->
## homoskedastic standard SE undercovers (robust/cluster/HAC needed) and the
## Bai/SWW B/N + C/T bias is enlarged.  Serial correlation -> standard & robust
## SE undercover; cluster (and HAC) SE correct it (SWW Case B).
## IMPORTANT: serial correlation is applied ONLY to the STATIC DGPs (1, 2),
## which have strictly exogenous regressors.  On the dynamic DGPs (3, 4) the
## regressor is x = y_{i,t-1}; serially correlated errors there would make the
## lag endogenous (E[x u] != 0), which no bias correction can fix.  So the
## dynamic DGPs stay serially uncorrelated regardless of MC_ERR (the serial
## flag is ignored for them), exactly as in Su, Wang & Wang (2025).
err_type <- Sys.getenv("MC_ERR", unset = "homo")
if (!err_type %in% c("homo", "hetero", "serial", "serialhet"))
  stop("MC_ERR must be one of homo/hetero/serial/serialhet (got '", err_type, "').")
HETERO <- err_type %in% c("hetero", "serialhet")   # heteroskedastic errors
SERIAL <- err_type %in% c("serial", "serialhet")   # MA(1) serial correlation (static only)
## SE type for bias-corrected estimator.  Must match the error structure:
##   homo       -> standard (valid and efficient)
##   hetero     -> robust (accounts for cross-sectional heteroskedasticity)
##   serial     -> cluster/hac (accounts for within-unit serial correlation)
##   serialhet  -> cluster/hac (accounts for both)
## Balanced DGPs (1, 3) use ife() which has standard/robust/cluster but no hac;
## cluster SE allows arbitrary within-unit correlation, so it covers serial+het.
## Unbalanced DGP 2 uses ife_unbalanced() which supports hac directly.
## DGP 4 always uses "hac" (hardcoded below, not controlled by BC_SE_*).
BC_SE_BAL <- if (SERIAL) "cluster" else if (HETERO) "robust" else "standard"
BC_SE_UNB <- if (SERIAL) "hac" else if (HETERO) "robust" else "standard"

## Reproducibility - unique seed per (dgp, N, TT, fill, b)
## Scheme: dgp*1e7 + N*1e4 + TT*1e2 + fill_int + b  (no collisions for realistic grids)
seed_base <- dgp_id * 10000000L + N * 10000L + TT * 100L + fill_int
RNGkind("L'Ecuyer-CMRG")

## Output location (override with env var XTIFE_MC_BASE if needed)
base_dir <- Sys.getenv("XTIFE_MC_BASE", unset = "/home/bc25911/xtife_mc")
out_dir  <- file.path(base_dir, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
## Each error structure gets a distinct filename tag so the variants never
## clobber each other and collate_mc.R can build separate table sets.
err_tag  <- switch(err_type, homo = "", hetero = "_het",
                   serial = "_ser", serialhet = "_serhet")
out_file <- file.path(out_dir,
  sprintf("dgp%d_N%d_T%d_f%d_%dreps%s.rds", dgp_id, N, TT, fill_int, B, err_tag))

cat(sprintf("[mc_xtife_hpc]  DGP=%d  N=%d  T=%d  fill=%.0f%%  cores=%d  B=%d  err=%s\n",
            dgp_id, N, TT, fill_int, nc, B, err_type))
cat(sprintf("  r_true=%d  r_max=%d  beta_true=%.1f  seed_base=%d\n",
            r_true, r_max, beta_true, seed_base))
cat(sprintf("  Output: %s\n", out_file))

# ==============================================================================
# 1 - DGP Generator Functions
# ==============================================================================

## ---------------------------------------------------------------------------
## DGP design notes (see data-raw/hpc/DECISIONS_mc_audit.md for the full audit).
##   * Loadings lambda_i ~ N(0, I_r), factors F_t ~ N(0, I_r): keeps the factor
##     structure well-conditioned so SVT / IC(BIC) factor selection is reliable.
##   * Errors u_it ~ N(0,1) i.i.d. homoskedastic.  Under homoskedasticity the
##     Bai (2009) / SWW2025 incidental-parameter bias is O(1/N+1/T) but small;
##     the bias correction is therefore a *small refinement* here (validated in
##     Part A to be correct and not to over-correct), not a dramatic effect.
##     Demonstrating a large correctable static bias requires factor-correlated
##     missingness (SWW Pattern 2) or strong heteroskedasticity, both of which
##     degrade factor selection / coverage in this pipeline (see DECISIONS doc)
##     and are deferred.
##   * Static regressor (DGP 1/2): strictly exogenous, x = gamma*(lambda'F)+eps.
##   * Dynamic regressor (DGP 3/4): x_it = y_{i,t-1} (predetermined), AR(1)
##     factors, SERIALLY UNCORRELATED errors (theta_ma = 0).  Serial errors with
##     a lagged-y regressor would be contemporaneously endogenous (E[x u] != 0),
##     which no bias correction can fix; this was the previous fatal MC bug.
## ---------------------------------------------------------------------------

## Heteroskedasticity design (used only when MC_ERR = "hetero").
## The error sd is tied to the common factor term fc = lambda_i'F_t:
##   sigma_it = 0.4 + 0.8 * |fc_it| / sd(fc),  rescaled so mean(sigma^2) = 1.
## This keeps the overall noise level comparable to the homoskedastic case but
## makes the error variance correlated with the regressor (which loads on the
## same factors) -> homoskedastic standard SE undercovers; robust/cluster/HAC
## SE are needed; and the Bai/SWW B/N + C/T bias is enlarged.  The computation
## is inlined in each generator (no helper function) so the heteroskedasticity
## logic always travels with the DGP and cannot go out of scope.

## DGP 1: static balanced
##   y_it = beta * x_it + lambda_i' F_t + u_it
##   x_it = (lambda_i + mu_i)' F_t + eps_it   [strictly exog, corr. with factors]
##   F_t, lambda_i, mu_i ~ N(0, I_r) i.i.d.;  eps_it ~ N(0,1)
##   u_it: homo -> N(0,1);  hetero -> sigma_it * N(0,1) with sigma_it tied to fc.
##   The extra loading mu_i (independent of lambda_i, as in SWW2025) makes x
##   correlated with y's factor term lambda_i'F_t (-> OLS bias, motivates IFE)
##   WITHOUT x being a clean combination of the factors -- otherwise x would
##   absorb one factor direction and IC(BIC) would select r = 1 instead of 2.
dgp_static_balanced <- function(N, TT, r, beta, het = HETERO, serial = SERIAL,
                                 rho_x = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  F_mat       <- matrix(rnorm(TT * r), TT, r)              # iid factors (all variants)
  Lambda_true <- matrix(rnorm(N  * r), N,  r)
  Mu          <- matrix(rnorm(N  * r), N,  r)        # independent extra loading
  fc          <- as.vector(Lambda_true %*% t(F_mat))       # lambda_i' F_t (y factor)
  xc          <- as.vector((Lambda_true + Mu) %*% t(F_mat))# (lambda+mu)' F_t
  df          <- expand.grid(unit = seq_len(N), time = seq_len(TT))
  ## Idiosyncratic regressor noise.  homo/hetero: i.i.d.  serial: AR(1) over time
  ## in the IDIOSYNCRATIC component (NOT the factor part).  The factor-driven part
  ## of x is removed by the IFE projection M_F, so to make serial correlation
  ## actually affect inference the serial part must live in the idiosyncratic
  ## component, which survives the projection.  (het/serial transforms below
  ## consume no RNG, so the homo path is byte-identical.)
  xn <- rnorm(N * TT)                                      # X noise (RNG position fixed)
  if (serial) {
    zm <- matrix(xn, nrow = N, ncol = TT)                 # zm[i,t], unit-fast
    for (tt in 2L:TT)
      zm[, tt] <- rho_x * zm[, tt - 1L] + sqrt(1 - rho_x^2) * zm[, tt]  # AR(1), var 1
    xn <- as.vector(zm)
  }
  df$X        <- xc + xn
  err         <- rnorm(N * TT)              # base innovations e_it (RNG position fixed)
  if (het) {                                # inline heteroskedasticity
    s   <- 0.4 + 0.8 * abs(fc) / stats::sd(fc)   # sigma_it tied to factor term
    err <- (s / sqrt(mean(s^2))) * err            # rescale to unit avg variance
  }
  if (serial) {                             # MA(1): v_it = (e_it + e_{i,t-1})/sqrt(2)
    em  <- matrix(err, nrow = N, ncol = TT)        # em[i,t] = error of unit i at time t
    em  <- (em + cbind(0, em[, -TT, drop = FALSE])) / sqrt(2)
    err <- as.vector(em)                           # back to unit-fast vector
  }
  df$Y        <- beta * df$X + fc + err
  df
}

## DGP 2: static unbalanced (Pattern 1 MCAR dropout from DGP 1)
dgp_static_unbalanced <- function(N, TT, r, beta, fill = 0.8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  df_full <- dgp_static_balanced(N, TT, r, beta, seed = NULL)
  keep    <- sort(sample(nrow(df_full), floor(fill * nrow(df_full))))
  df_full[keep, ]
}

## DGP 3: dynamic balanced (Moon & Weidner 2017)
##   y_it = beta * y_{i,t-1} + lambda_i' F_t + u_it
##   F_t  = rho_f * F_{t-1} + eta_t  (AR(1), stationary unit variance)
##   u_it ~ N(0,1) i.i.d. (serially uncorrelated -> x = lag(y) predetermined)
dgp_dynamic_balanced <- function(N, TT, r, beta,
                                  rho_f = 0.5, sigma_u = 1.0,
                                  T0_burnin = 50L, het = HETERO, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  TT_total <- TT + T0_burnin
  sv_eta   <- sqrt(1 - rho_f^2)
  F_mat    <- matrix(0, TT_total, r)
  F_mat[1L, ] <- rnorm(r)
  for (tt in 2L:TT_total)
    F_mat[tt, ] <- rho_f * F_mat[tt - 1L, ] + rnorm(r, sd = sv_eta)
  Lambda_true <- matrix(rnorm(N * r), N, r)
  ## NOTE: errors are kept serially UNCORRELATED here regardless of MC_ERR
  ## ("serial"/"serialhet" do not add MA to the dynamic DGPs): with x = y_{i,t-1}
  ## a serially correlated error would be contemporaneously endogenous.
  u_mat  <- matrix(rnorm(TT_total * N), TT_total, N)   # serially uncorrelated
  if (het) {                                           # tie sigma to factor term
    lf_all <- F_mat %*% t(Lambda_true)                 # TT_total x N : lambda_i'F_t
    sig    <- 0.4 + 0.8 * abs(lf_all) / stats::sd(as.vector(lf_all))
    u_mat  <- (sig / sqrt(mean(sig^2))) * u_mat
  }
  Y_mat  <- matrix(0, TT_total, N)
  lf1    <- as.vector(F_mat[1L, , drop = FALSE] %*% t(Lambda_true))
  Y_mat[1L, ] <- lf1 + sigma_u * u_mat[1L, ]
  for (tt in 2L:TT_total) {
    ft <- as.vector(F_mat[tt, , drop = FALSE] %*% t(Lambda_true))
    Y_mat[tt, ] <- beta * Y_mat[tt - 1L, ] + ft + sigma_u * u_mat[tt, ]
  }
  idx_y <- (T0_burnin + 2L):TT_total
  idx_x <- (T0_burnin + 1L):(TT_total - 1L)
  n_t   <- length(idx_y)
  df <- data.frame(
    unit = rep(seq_len(N),  each  = n_t),
    time = rep(seq_len(n_t), times = N),
    Y    = as.vector(Y_mat[idx_y, ]),
    X    = as.vector(Y_mat[idx_x, ])
  )
  list(df = df, TT_eff = n_t)
}

## DGP 4: dynamic unbalanced (Pattern 1 MCAR dropout from DGP 3)
dgp_dynamic_unbalanced <- function(N, TT, r, beta,
                                    fill = 0.8, ..., seed = NULL) {
  res     <- dgp_dynamic_balanced(N, TT, r, beta, ..., seed = seed)
  df_full <- res$df
  keep    <- sort(sample(nrow(df_full), floor(fill * nrow(df_full))))
  list(df = df_full[keep, ], TT_eff = res$TT_eff)
}

# ==============================================================================
# 2 - One-Replication Estimation
# ==============================================================================

## Column names of the output vector (21 statistics per rep)
OUT_NAMES <- c(
  "b_ols",                                    # OLS (r=0)
  "b_kn",  "b_bc",                            # IFE r-known: raw, bias-corrected
  "se_kn_std", "se_kn_rob", "se_kn_cl",      # standard errors (r known)
  "cov_kn_std", "cov_kn_rob", "cov_kn_cl",   # 95% CI coverage (r known)
  "cov_bc_std",                               # coverage of BC estimator (std SE)
  "size_kn_std",                              # size: |t_std| > 1.96 under H0
  "r_hat", "r_correct",                       # factor selection
  "b_un",  "cov_un_std",                      # r-unknown two-step estimator
  "conv",  "n_iter",                          # convergence info
  "b_hac", "cov_kn_hac", "se_kn_hac",        # HAC statistics (DGPs 3 & 4)
  "conv_sel"                                  # factor selection converged
)
N_OUT <- length(OUT_NAMES)   # = 21

run_one_rep <- function(N, TT, r_true, r_max, dgp_id,
                        beta_true, fill = 0.8, seed = NULL) {
  ## Generate data
  df <- switch(as.character(dgp_id),
    "1" = dgp_static_balanced(N, TT, r_true, beta_true, seed = seed),
    "2" = dgp_static_unbalanced(N, TT, r_true, beta_true,
                                 fill = fill, seed = seed),
    "3" = dgp_dynamic_balanced(N, TT, r_true, beta_true, seed = seed)$df,
    "4" = dgp_dynamic_unbalanced(N, TT, r_true, beta_true,
                                  fill = fill, seed = seed)$df
  )
  idx <- c("unit", "time")
  out <- rep(NA_real_, N_OUT)
  names(out) <- OUT_NAMES
  .err_msg <- NA_character_   # captured on crash so the master can surface it

  tryCatch({

    ## (A) OLS baseline
    ## Note: force="none" throughout because all four DGPs have NO additive
    ## unit/time fixed effects — the only confounders are the interactive factors
    ## lambda_i'f_t.  For DGPs 1/2 this creates the classic Bai (2009) factor-
    ## omission bias.  For DGPs 3/4 (dynamic), this is OLS ignoring BOTH the
    ## factor structure AND the lagged-Y endogeneity — the largest possible bias,
    ## demonstrating why IFE is needed.  The canonical Nickell (1981) benchmark
    ## (force="unit") is not applicable here since there are no unit FEs in the DGP.
    ## OLS baselines are specified identically across all four DGPs:
    ## ife(r=0, force="none") grand-mean-demeans Y and X (Y - mu_Y), which by
    ## Frisch-Waugh-Lovell is exactly pooled OLS WITH an intercept. The lm()
    ## branch for the unbalanced DGPs therefore also keeps the intercept
    ## (lm(Y ~ X)) so that all four OLS columns measure the same estimand.
    out["b_ols"] <- if (dgp_id %in% c(1L, 3L)) {
      fit0 <- ife(Y ~ X, data = df, index = idx,
                  r = 0L, force = "none", se = "standard")
      unname(fit0$coef["X"])
    } else {
      unname(coef(lm(Y ~ X, data = df))["X"])
    }

    ## (B) IFE r-known
    if (dgp_id == 1L) {
      fit_std <- ife(Y ~ X, data = df, index = idx, r = r_true,
                     force = "none", se = "standard", bias_corr = FALSE)
      fit_bc  <- ife(Y ~ X, data = df, index = idx, r = r_true,
                     force = "none", se = BC_SE_BAL, bias_corr = TRUE)
      fit_rob <- ife(Y ~ X, data = df, index = idx, r = r_true,
                     force = "none", se = "robust",   bias_corr = FALSE)
      fit_cl  <- ife(Y ~ X, data = df, index = idx, r = r_true,
                     force = "none", se = "cluster",  bias_corr = FALSE)

    } else if (dgp_id == 2L) {
      fit_std <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "standard", bias_corr = FALSE, init = "nnr"))
      fit_bc  <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = BC_SE_UNB, bias_corr = TRUE, exog = "strict", init = "nnr"))
      fit_rob <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "robust",   bias_corr = FALSE, init = "nnr"))
      fit_cl  <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "cluster",  bias_corr = FALSE, init = "nnr"))

    } else if (dgp_id == 3L) {
      fit_std <- ife(Y ~ X, data = df, index = idx, r = r_true, force = "none",
                     se = "standard", method = "dynamic", bias_corr = FALSE)
      fit_bc  <- ife(Y ~ X, data = df, index = idx, r = r_true, force = "none",
                     se = BC_SE_BAL, method = "dynamic", bias_corr = TRUE, M1 = 1L)
      fit_rob <- ife(Y ~ X, data = df, index = idx, r = r_true, force = "none",
                     se = "robust",   method = "dynamic", bias_corr = FALSE)
      fit_cl  <- ife(Y ~ X, data = df, index = idx, r = r_true, force = "none",
                     se = "cluster",  method = "dynamic", bias_corr = FALSE)

    } else {   # dgp_id == 4L
      fit_std <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "standard", bias_corr = FALSE, exog = "strict", init = "nnr"))
      fit_bc  <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "hac", bias_corr = TRUE, exog = "weak", init = "nnr"))
      fit_rob <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "robust",   bias_corr = FALSE, exog = "strict", init = "nnr"))
      fit_cl  <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "cluster",  bias_corr = FALSE, exog = "strict", init = "nnr"))
    }

    b_kn  <- unname(fit_std$coef["X"])
    b_bc  <- unname(fit_bc$coef["X"])
    se_std <- unname(fit_std$se["X"])
    se_rob <- unname(fit_rob$se["X"])
    se_cl  <- unname(fit_cl$se["X"])
    ci_std <- fit_std$ci["X", ];  ci_rob <- fit_rob$ci["X", ]
    ci_cl  <- fit_cl$ci["X",  ];  ci_bc  <- fit_bc$ci["X",  ]
    conv   <- as.integer(fit_std$converged)
    n_iter <- if (!is.null(fit_std$n_iter)) fit_std$n_iter else NA_integer_

    ## HAC / cluster coverage
    se_hac <- NA_real_; ci_hac <- c(NA_real_, NA_real_)
    if (dgp_id == 2L) {
      ## Static unbalanced: HAC SE (Bartlett kernel) on the raw IFE, strictly
      ## exogenous (SWW Case B).  Under serial correlation the homoskedastic
      ## standard SE undercovers; this HAC column shows the correction.  Under
      ## homo/hetero (no serial correlation) HAC ~ robust, so it is harmless.
      fit_hac_raw <- tryCatch(suppressWarnings(
        ife_unbalanced(Y ~ X, data = df, index = idx,
                       r = r_true, se = "hac", bias_corr = FALSE,
                       exog = "strict", init = "nnr")),
        error = function(e) NULL)
      if (!is.null(fit_hac_raw)) {
        se_hac <- unname(fit_hac_raw$se["X"])
        ci_hac <- fit_hac_raw$ci["X", ]
      }

    } else if (dgp_id == 3L) {
      ## Balanced dynamic: cluster SE is the best available proxy for HAC
      ## (ife() has no HAC option; cluster accounts for all within-unit correlation)
      se_hac <- se_cl; ci_hac <- ci_cl

    } else if (dgp_id == 4L) {
      ## Raw IFE with HAC SE (NO bias correction) — strictly exog SE formula applied
      ## to a weakly exog DGP.  This shows HAC SE alone is insufficient: it handles
      ## serial correlation but NOT the dynamic (weak exog) bias.
      ## Contrast with BC.Cov = fit_bc (BC + HAC), which should give ~95%.
      fit_hac_raw <- tryCatch(suppressWarnings(
        ife_unbalanced(Y ~ X, data = df, index = idx,
                       r = r_true, se = "hac", bias_corr = FALSE,
                       exog = "strict", init = "nnr")),
        error = function(e) NULL)
      if (!is.null(fit_hac_raw)) {
        se_hac <- unname(fit_hac_raw$se["X"])
        ci_hac <- fit_hac_raw$ci["X", ]
      }
    }

    out["b_kn"]       <- b_kn
    out["b_bc"]       <- b_bc
    out["se_kn_std"]  <- se_std
    out["se_kn_rob"]  <- se_rob
    out["se_kn_cl"]   <- se_cl
    out["se_kn_hac"]  <- se_hac
    out["b_hac"]      <- b_kn   # always raw IFE (b_bc tracked separately)
    out["cov_kn_std"] <- as.integer(ci_std[1L] <= beta_true & beta_true <= ci_std[2L])
    out["cov_kn_rob"] <- as.integer(ci_rob[1L] <= beta_true & beta_true <= ci_rob[2L])
    out["cov_kn_cl"]  <- as.integer(ci_cl[1L]  <= beta_true & beta_true <= ci_cl[2L])
    ## Note: cov_bc_std stores the coverage of the BC estimator with the SE type
    ## that matches the error structure.  Balanced DGPs 1/3 use se = BC_SE_BAL:
    ##   homo -> "standard",  hetero -> "robust",  serial/serialhet -> "cluster".
    ## Unbalanced DGP 2 uses se = BC_SE_UNB:
    ##   homo -> "standard",  hetero/serial/serialhet -> "hac".
    ## DGP 4 always uses se = "hac" (dynamic unbalanced needs HAC regardless).
    ## collate_mc.R relabels the column appropriately per error type.
    out["cov_bc_std"] <- as.integer(ci_bc[1L]  <= beta_true & beta_true <= ci_bc[2L])
    out["cov_kn_hac"] <- if (!is.na(ci_hac[1L]))
                           as.integer(ci_hac[1L] <= beta_true & beta_true <= ci_hac[2L])
                         else NA_real_
    out["size_kn_std"] <- as.integer(abs((b_kn - beta_true) / se_std) > 1.96)
    out["conv"]   <- conv
    out["n_iter"] <- n_iter

    ## (C) Factor selection (r unknown)
    r_hat <- tryCatch({
      if (dgp_id == 1L || dgp_id == 3L) {
        sel <- suppressWarnings(
          ife_select_r(Y ~ X, data = df, index = idx,
                       r_max = r_max, force = "none", verbose = FALSE))
        as.integer(attr(sel, "suggested")["IC_bic"])
      } else {
        sel <- suppressWarnings(
          ife_select_r_unb(Y ~ X, data = df, index = idx, verbose = FALSE))
        as.integer(sel$r_hat)
      }
    }, error = function(e) NA_integer_)

    out["conv_sel"]  <- as.integer(!is.na(r_hat))
    r_hat            <- if (is.na(r_hat)) r_true else max(1L, min(r_hat, r_max))
    out["r_hat"]     <- r_hat
    out["r_correct"] <- as.integer(r_hat == r_true)

    ## (D) Two-step estimator (r unknown, then BC)
    fit_un <- tryCatch(suppressWarnings({
      if (dgp_id == 1L)
        ife(Y ~ X, data = df, index = idx, r = r_hat, force = "none",
            se = "standard", bias_corr = TRUE)
      else if (dgp_id == 2L)
        ife_unbalanced(Y ~ X, data = df, index = idx, r = r_hat,
                       se = "standard", bias_corr = TRUE, exog = "strict", init = "nnr")
      else if (dgp_id == 3L)
        ife(Y ~ X, data = df, index = idx, r = r_hat, force = "none",
            se = "standard", method = "dynamic", bias_corr = TRUE, M1 = 1L)
      else
        ife_unbalanced(Y ~ X, data = df, index = idx, r = r_hat,
                       se = "hac", bias_corr = TRUE, exog = "weak", init = "nnr")
    }), error = function(e) NULL)

    if (!is.null(fit_un)) {
      ci_un <- fit_un$ci["X", ]
      out["b_un"]      <- unname(fit_un$coef["X"])
      ## Note: for DGP 4, fit_un uses se="hac", so ci_un is a HAC CI.
      ## cov_un_std therefore stores two-step BC+HAC coverage for DGP 4.
      out["cov_un_std"] <- as.integer(ci_un[1L] <= beta_true & beta_true <= ci_un[2L])
    }

  }, error = function(e) .err_msg <<- conditionMessage(e))   # capture, leave NAs

  if (!is.na(.err_msg)) attr(out, "err") <- .err_msg   # survives back to master
  out
}

# ==============================================================================
# 3 - Summary Statistics
# ==============================================================================

summarise_reps <- function(mat, beta_true, r_true) {
  cov_mean <- function(x) mean(x, na.rm = TRUE)
  cov_se   <- function(x) {
    m <- cov_mean(x)
    sqrt(m * (1 - m) / max(1L, sum(!is.na(x))))
  }
  b_ols <- mat[, "b_ols"];  b_kn <- mat[, "b_kn"]
  b_bc  <- mat[, "b_bc"];   b_un <- mat[, "b_un"]
  se_std <- mat[, "se_kn_std"]

  data.frame(
    OLS.Bias    = mean(b_ols - beta_true, na.rm = TRUE),
    IFE.Bias    = mean(b_kn  - beta_true, na.rm = TRUE),
    IFE.SD      = sd(b_kn,   na.rm = TRUE),
    IFE.RMSE    = sqrt(mean((b_kn - beta_true)^2, na.rm = TRUE)),
    IFE.SE_SD   = mean(se_std, na.rm = TRUE) / sd(b_kn, na.rm = TRUE),
    IFE.Cov.Std = cov_mean(mat[, "cov_kn_std"]),
    IFE.Cov.Rob = cov_mean(mat[, "cov_kn_rob"]),
    IFE.Cov.Cl  = cov_mean(mat[, "cov_kn_cl"]),
    IFE.Cov.HAC = cov_mean(mat[, "cov_kn_hac"]),
    IFE.Cov.SE  = cov_se(mat[, "cov_kn_std"]),
    IFE.Size    = mean(mat[, "size_kn_std"], na.rm = TRUE),
    BC.Bias     = mean(b_bc - beta_true, na.rm = TRUE),
    BC.SD       = sd(b_bc,  na.rm = TRUE),
    BC.RMSE     = sqrt(mean((b_bc - beta_true)^2, na.rm = TRUE)),
    BC.Cov.Std  = cov_mean(mat[, "cov_bc_std"]),
    BC.Cov.SE   = cov_se(mat[, "cov_bc_std"]),
    r.mean      = mean(mat[, "r_hat"],     na.rm = TRUE),
    r.correct   = mean(mat[, "r_correct"], na.rm = TRUE),
    r.under     = mean(mat[, "r_hat"] < r_true,  na.rm = TRUE),
    r.over      = mean(mat[, "r_hat"] > r_true,  na.rm = TRUE),
    UN.Bias     = mean(b_un - beta_true, na.rm = TRUE),
    UN.SD       = sd(b_un,  na.rm = TRUE),
    UN.RMSE     = sqrt(mean((b_un - beta_true)^2, na.rm = TRUE)),
    UN.Cov.Std  = cov_mean(mat[, "cov_un_std"]),
    conv.rate   = mean(mat[, "conv"],    na.rm = TRUE),
    n.iter.med  = median(mat[, "n_iter"], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# 4 - Run simulation for this cell
# ==============================================================================

## Disable multi-threaded BLAS inside forked workers (deadlock prevention)
Sys.setenv(OPENBLAS_NUM_THREADS   = "1",
           OMP_NUM_THREADS        = "1",
           MKL_NUM_THREADS        = "1",
           VECLIB_MAXIMUM_THREADS = "1")

seeds <- seed_base + seq_len(B)

cat(sprintf("[%s]  Dispatching %d replications on %d cores ...\n",
            format(Sys.time(), "%H:%M:%S"), B, nc))
t_start <- proc.time()["elapsed"]

reps <- tryCatch(
  parallel::mclapply(
    seq_len(B),
    function(b) run_one_rep(N, TT, r_true, r_max, dgp_id,
                             beta_true, fill, seed = seeds[b]),
    mc.cores    = nc,
    mc.set.seed = FALSE
  ),
  error = function(e) {
    cat("  [WARNING] mclapply failed - falling back to sequential lapply\n")
    lapply(seq_len(B),
           function(b) run_one_rep(N, TT, r_true, r_max, dgp_id,
                                    beta_true, fill, seed = seeds[b]))
  }
)

## Surface per-rep error messages (captured as attr "err") BEFORE rbind drops
## attributes.  If reps crashed, print the most common error so failures are
## diagnosable from the .out log instead of showing silent NaNs.
rep_errs <- unlist(lapply(reps, function(x) attr(x, "err")), use.names = FALSE)
if (length(rep_errs) > 0L) {
  tab <- sort(table(rep_errs), decreasing = TRUE)
  cat(sprintf("  [DIAG] %d/%d reps errored.  Most common message(s):\n", length(rep_errs), B))
  for (j in seq_len(min(3L, length(tab))))
    cat(sprintf("    [%d x] %s\n", tab[j], names(tab)[j]))
}

## Replace crashed worker output (try-error / NULL / wrong length) with NAs
na_vec <- { v <- rep(NA_real_, N_OUT); names(v) <- OUT_NAMES; v }
reps <- lapply(reps, function(x) {
  if (inherits(x, "try-error") || is.null(x) || length(x) != N_OUT) na_vec else x
})

raw_mat <- do.call(rbind, reps)
summ    <- summarise_reps(raw_mat, beta_true, r_true)
summ$N  <- N;  summ$TT <- TT;  summ$fill <- fill;  summ$dgp_id <- dgp_id

t_elapsed <- proc.time()["elapsed"] - t_start
n_conv    <- sum(raw_mat[, "conv"] == 1L, na.rm = TRUE)
n_na      <- sum(is.na(raw_mat[, "b_kn"]))

cat(sprintf("[%s]  Done in %.1f s  (%.0f ms/rep)\n",
            format(Sys.time(), "%H:%M:%S"), t_elapsed, 1000 * t_elapsed / B))
cat(sprintf("  Converged: %d/%d   NA reps: %d\n", n_conv, B, n_na))
cat(sprintf("  IFE bias = %.4f   BC bias = %.4f\n",
            summ$IFE.Bias, summ$BC.Bias))
cat(sprintf("  Coverage: raw std = %.3f  raw rob = %.3f  raw cl = %.3f  |  BC = %.3f\n",
            summ$IFE.Cov.Std, summ$IFE.Cov.Rob, summ$IFE.Cov.Cl, summ$BC.Cov.Std))

# ==============================================================================
# 5 - Save results
# ==============================================================================

saveRDS(
  list(
    params  = list(dgp_id    = dgp_id,
                   N         = N,
                   TT        = TT,
                   fill      = fill,
                   B         = B,
                   r_true    = r_true,
                   r_max     = r_max,
                   beta_true = beta_true,
                   nc        = nc,
                   seed_base = seed_base),
    summ    = summ,
    raw_mat = raw_mat,      # B x 21 matrix (168 KB for B=1000; kept for diagnostics)
    timing  = list(elapsed_sec  = t_elapsed,
                   ms_per_rep   = 1000 * t_elapsed / B,
                   n_converged  = n_conv,
                   n_na_reps    = n_na)
  ),
  file     = out_file,
  compress = "gzip"
)

cat(sprintf("  Saved: %s\n", out_file))
