# ==============================================================================
# mc_xtife_hpc.R
# HPC worker script: runs ONE (dgp_id, N, TT, fill) parameter cell
# and saves results to data-raw/hpc/results/.
#
# Usage (called by run_mc_job.sh):
#   Rscript mc_xtife_hpc.R <dgp_id> <N> <TT> <fill×100> <ncores> [B]
#
#   dgp_id  : 1 = static balanced,  2 = static unbalanced,
#             3 = dynamic balanced,  4 = dynamic unbalanced
#   N       : number of cross-sectional units
#   TT      : number of time periods
#   fill×100: e.g. 80 means 80 % of cells observed (used for DGPs 2 & 4)
#   ncores  : parallel workers ($NSLOTS from SGE)
#   B       : replications (default 1000)
#
# Output: $XTIFE_MC_BASE/results/dgp{d}_N{N}_T{TT}_f{fill}_{B}reps.rds
#   List: params, summ (summary stats), raw_mat (B×21 matrix), timing
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
    "Option A — install to user library on HPC:\n",
    "  mkdir -p ~/R/library\n",
    "  R -e 'install.packages(\"xtife\", lib=\"~/R/library\")'\n",
    "  Then add: export R_LIBS_USER=~/R/library to run_mc_job.sh\n\n",
    "Option B — source from files (no installation needed):\n",
    "  rsync -av R/ife.R R/ife_unbalanced.R ",
    "bc25911@hpc.essex.ac.uk:", .xtife_src, "/\n"
  )
}
rm(.xtife_base, .xtife_src, list = ls(pattern = "^\\.f$"))

# ==============================================================================
# 0 — Command-line arguments & fixed parameters
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L)
  stop("Usage: Rscript mc_xtife_hpc.R <dgp_id> <N> <TT> <fill×100> <ncores> [B]")

dgp_id   <- as.integer(args[1L])
N        <- as.integer(args[2L])
TT       <- as.integer(args[3L])
fill_int <- as.integer(args[4L])    # e.g. 80
nc       <- as.integer(args[5L])    # $NSLOTS
B        <- if (length(args) >= 6L) as.integer(args[6L]) else 1000L
fill     <- fill_int / 100          # e.g. 0.80

## Fixed across all cells
r_true    <- 2L
r_max     <- 5L
beta_true <- 1.0

## Reproducibility — unique seed per (dgp, N, TT, fill, b)
## Scheme: dgp*1e7 + N*1e4 + TT*1e2 + fill_int + b  (no collisions for realistic grids)
seed_base <- dgp_id * 10000000L + N * 10000L + TT * 100L + fill_int
RNGkind("L'Ecuyer-CMRG")

## Output location (override with env var XTIFE_MC_BASE if needed)
base_dir <- Sys.getenv("XTIFE_MC_BASE", unset = "/home/bc25911/xtife_mc")
out_dir  <- file.path(base_dir, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir,
  sprintf("dgp%d_N%d_T%d_f%d_%dreps.rds", dgp_id, N, TT, fill_int, B))

cat(sprintf("[mc_xtife_hpc]  DGP=%d  N=%d  T=%d  fill=%.0f%%  cores=%d  B=%d\n",
            dgp_id, N, TT, fill_int, nc, B))
cat(sprintf("  r_true=%d  r_max=%d  beta_true=%.1f  seed_base=%d\n",
            r_true, r_max, beta_true, seed_base))
cat(sprintf("  Output: %s\n", out_file))

# ==============================================================================
# 1 — DGP Generator Functions
# ==============================================================================

## DGP 1: static balanced
##   y_it = beta * x_it + lambda_i' F_t + u_it
##   x_it = gamma * (lambda_i' F_t) + eps_it   (endogeneity -> TWFE biased)
dgp_static_balanced <- function(N, TT, r, beta,
                                 gamma = 1.0, sigma_u = 1.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  F_true      <- matrix(rnorm(TT * r), TT, r)
  Lambda_true <- matrix(rnorm(N  * r), N,  r)
  df          <- expand.grid(unit = seq_len(N), time = seq_len(TT))
  fc          <- as.vector(Lambda_true %*% t(F_true))   # N*TT, unit-fast col-major
  df$X <- gamma * fc + rnorm(N * TT)
  df$Y <- beta  * df$X + fc + rnorm(N * TT, sd = sigma_u)
  df
}

## DGP 2: static unbalanced (MCAR dropout from DGP 1)
dgp_static_unbalanced <- function(N, TT, r, beta,
                                   gamma = 1.0, sigma_u = 1.0,
                                   fill = 0.8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  df_full <- dgp_static_balanced(N, TT, r, beta, gamma, sigma_u, seed = NULL)
  keep    <- sort(sample(nrow(df_full), floor(fill * nrow(df_full))))
  df_full[keep, ]
}

## DGP 3: dynamic balanced (Moon & Weidner 2017)
##   y_it = beta * y_{i,t-1} + lambda_i' F_t + u_it
##   F_t  = rho_f * F_{t-1} + eta_t  (AR(1), unit variance)
##   u_it = (e_it + theta_ma * e_{i,t-1}) / sqrt(1+theta^2)  (MA(1))
dgp_dynamic_balanced <- function(N, TT, r, beta,
                                  rho_f = 0.5, theta_ma = 0.5, sigma_u = 1.0,
                                  T0_burnin = 50L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  TT_total <- TT + T0_burnin
  sv_eta   <- sqrt(1 - rho_f^2)
  F_mat    <- matrix(0, TT_total, r)
  F_mat[1L, ] <- rnorm(r)
  for (tt in 2L:TT_total)
    F_mat[tt, ] <- rho_f * F_mat[tt - 1L, ] + rnorm(r, sd = sv_eta)
  Lambda_true <- matrix(rnorm(N * r), N, r)
  e_mat  <- matrix(rnorm(TT_total * N), TT_total, N)
  sv_ma  <- sqrt(1 + theta_ma^2)
  u_mat  <- (e_mat + theta_ma * rbind(matrix(0, 1L, N), e_mat[-TT_total, ])) / sv_ma
  Y_mat  <- matrix(0, TT_total, N)
  Y_mat[1L, ] <- as.vector(F_mat[1L, , drop = FALSE] %*% t(Lambda_true)) +
                 sigma_u * u_mat[1L, ]
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

## DGP 4: dynamic unbalanced (MCAR dropout from DGP 3)
dgp_dynamic_unbalanced <- function(N, TT, r, beta,
                                    fill = 0.8, ..., seed = NULL) {
  res     <- dgp_dynamic_balanced(N, TT, r, beta, ..., seed = seed)
  df_full <- res$df
  keep    <- sort(sample(nrow(df_full), floor(fill * nrow(df_full))))
  list(df = df_full[keep, ], TT_eff = res$TT_eff)
}

# ==============================================================================
# 2 — One-Replication Estimation
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

  tryCatch({

    ## (A) OLS baseline
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
                     force = "none", se = "standard", bias_corr = TRUE)
      fit_rob <- ife(Y ~ X, data = df, index = idx, r = r_true,
                     force = "none", se = "robust",   bias_corr = FALSE)
      fit_cl  <- ife(Y ~ X, data = df, index = idx, r = r_true,
                     force = "none", se = "cluster",  bias_corr = FALSE)

    } else if (dgp_id == 2L) {
      fit_std <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "standard", bias_corr = FALSE))
      fit_bc  <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "standard", bias_corr = TRUE, exog = "strict"))
      fit_rob <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "robust",   bias_corr = FALSE))
      fit_cl  <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "cluster",  bias_corr = FALSE))

    } else if (dgp_id == 3L) {
      fit_std <- ife(Y ~ X, data = df, index = idx, r = r_true, force = "none",
                     se = "standard", method = "dynamic", bias_corr = FALSE)
      fit_bc  <- ife(Y ~ X, data = df, index = idx, r = r_true, force = "none",
                     se = "standard", method = "dynamic", bias_corr = TRUE, M1 = 1L)
      fit_rob <- ife(Y ~ X, data = df, index = idx, r = r_true, force = "none",
                     se = "robust",   method = "dynamic", bias_corr = FALSE)
      fit_cl  <- ife(Y ~ X, data = df, index = idx, r = r_true, force = "none",
                     se = "cluster",  method = "dynamic", bias_corr = FALSE)

    } else {   # dgp_id == 4L
      fit_std <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "standard", bias_corr = FALSE, exog = "strict"))
      fit_bc  <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "hac", bias_corr = TRUE, exog = "weak"))
      fit_rob <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "robust",   bias_corr = FALSE, exog = "strict"))
      fit_cl  <- suppressWarnings(ife_unbalanced(Y ~ X, data = df, index = idx,
                   r = r_true, se = "cluster",  bias_corr = FALSE, exog = "strict"))
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

    ## HAC / cluster coverage (DGPs 3 and 4)
    se_hac <- NA_real_; ci_hac <- c(NA_real_, NA_real_)
    if (dgp_id == 3L) {
      se_hac <- se_cl; ci_hac <- ci_cl      # cluster as proxy for dynamic balanced
    } else if (dgp_id == 4L) {
      se_hac <- unname(fit_bc$se["X"]); ci_hac <- ci_bc
    }

    out["b_kn"]       <- b_kn
    out["b_bc"]       <- b_bc
    out["se_kn_std"]  <- se_std
    out["se_kn_rob"]  <- se_rob
    out["se_kn_cl"]   <- se_cl
    out["se_kn_hac"]  <- se_hac
    out["b_hac"]      <- if (dgp_id == 4L) b_bc else b_kn
    out["cov_kn_std"] <- as.integer(ci_std[1L] <= beta_true & beta_true <= ci_std[2L])
    out["cov_kn_rob"] <- as.integer(ci_rob[1L] <= beta_true & beta_true <= ci_rob[2L])
    out["cov_kn_cl"]  <- as.integer(ci_cl[1L]  <= beta_true & beta_true <= ci_cl[2L])
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
                       se = "standard", bias_corr = TRUE, exog = "strict")
      else if (dgp_id == 3L)
        ife(Y ~ X, data = df, index = idx, r = r_hat, force = "none",
            se = "standard", method = "dynamic", bias_corr = TRUE, M1 = 1L)
      else
        ife_unbalanced(Y ~ X, data = df, index = idx, r = r_hat,
                       se = "hac", bias_corr = TRUE, exog = "weak")
    }), error = function(e) NULL)

    if (!is.null(fit_un)) {
      ci_un <- fit_un$ci["X", ]
      out["b_un"]      <- unname(fit_un$coef["X"])
      out["cov_un_std"] <- as.integer(ci_un[1L] <= beta_true & beta_true <= ci_un[2L])
    }

  }, error = function(e) invisible(NULL))   # leave all entries NA on crash

  out
}

# ==============================================================================
# 3 — Summary Statistics
# ==============================================================================

summarise_reps <- function(mat, beta_true) {
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
# 4 — Run simulation for this cell
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
    cat("  [WARNING] mclapply failed — falling back to sequential lapply\n")
    lapply(seq_len(B),
           function(b) run_one_rep(N, TT, r_true, r_max, dgp_id,
                                    beta_true, fill, seed = seeds[b]))
  }
)

## Replace crashed worker output (try-error / NULL / wrong length) with NAs
na_vec <- { v <- rep(NA_real_, N_OUT); names(v) <- OUT_NAMES; v }
reps <- lapply(reps, function(x) {
  if (inherits(x, "try-error") || is.null(x) || length(x) != N_OUT) na_vec else x
})

raw_mat <- do.call(rbind, reps)
summ    <- summarise_reps(raw_mat, beta_true)
summ$N  <- N;  summ$TT <- TT;  summ$fill <- fill;  summ$dgp_id <- dgp_id

t_elapsed <- proc.time()["elapsed"] - t_start
n_conv    <- sum(raw_mat[, "conv"] == 1L, na.rm = TRUE)
n_na      <- sum(is.na(raw_mat[, "b_kn"]))

cat(sprintf("[%s]  Done in %.1f s  (%.0f ms/rep)\n",
            format(Sys.time(), "%H:%M:%S"), t_elapsed, 1000 * t_elapsed / B))
cat(sprintf("  Converged: %d/%d   NA reps: %d\n", n_conv, B, n_na))
cat(sprintf("  IFE bias = %.4f   BC bias = %.4f   Coverage(std) = %.3f\n",
            summ$IFE.Bias, summ$BC.Bias, summ$IFE.Cov.Std))

# ==============================================================================
# 5 — Save results
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
    raw_mat = raw_mat,      # B × 21 matrix (168 KB for B=1000; kept for diagnostics)
    timing  = list(elapsed_sec  = t_elapsed,
                   ms_per_rep   = 1000 * t_elapsed / B,
                   n_converged  = n_conv,
                   n_na_reps    = n_na)
  ),
  file     = out_file,
  compress = "gzip"
)

cat(sprintf("  Saved: %s\n", out_file))
