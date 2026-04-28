# ============================================================================
# ife_unbalanced.R — Unbalanced Panel IFE via EM Algorithm
#
# Model: Y_it = X_it' beta + lambda_i' F_t + u_it
# No additive fixed effects (force = "none" only).
# Algorithm: Bai (2009) Appendix B EM procedure.
#
# Functions (internal → public):
#   .ife_em_inner()    — inner EM loop: estimates F and Lambda given beta
#   .ife_fit_unb()     — outer loop: alternates beta ↔ (F, Lambda)
#   .ife_se_unb()      — sandwich SE for unbalanced panel
#   ife_unbalanced()   — public wrapper, returns class "ife_unb"
#   print.ife_unb()    — print method
# ============================================================================


# ----------------------------------------------------------------------------
# .ife_em_inner — Inner EM loop (Bai 2009 Appendix B)
#
# Given a fixed beta (already absorbed into W_mat), iterates between:
#   E-step: impute missing cells with current lambda_i' F_t
#   M-step: SVD/eigen of the completed matrix to update F and Lambda
#
# @param W_mat       T x N matrix: observed W_it = Y_it - X_it'beta for
#                    observed cells; 0 elsewhere (caller sets this)
# @param obs_mask    T x N logical matrix: TRUE = observed cell
# @param Lambda_prev N x r loading matrix from previous outer iteration
#                    (used to warm-start imputation; zeros on first call)
# @param F_prev      T x r factor matrix from previous outer iteration
# @param N           number of units
# @param TT          number of time periods (= max T_i)
# @param r           number of factors
# @param tol_em      convergence tolerance on max fitted-value change at
#                    observed cells (default 1e-7)
# @param max_iter_em maximum EM inner iterations (default 500L)
#
# @return list(F_hat, Lambda_hat, n_iter_em, converged_em)
# ----------------------------------------------------------------------------
.ife_em_inner <- function(W_mat, obs_mask, Lambda_prev, F_prev,
                          N, TT, r,
                          tol_em      = 1e-7,
                          max_iter_em = 500L) {

  # Working completed matrix: start with observed values
  W_fill <- W_mat   # T x N, missing cells already 0

  # E-step initialisation: impute missing cells with previous fitted values
  # fitted_prev[t, i] = lambda_i^(h-1)' F_t^(h-1)
  if (any(r > 0L)) {
    fitted_prev <- F_prev %*% t(Lambda_prev)   # T x N
    W_fill[!obs_mask] <- fitted_prev[!obs_mask]
  }

  F_new      <- F_prev
  Lambda_new <- Lambda_prev

  for (h in seq_len(max_iter_em)) {

    # ---- M-step: update F via eigen of W_fill W_fill' / (N*TT) ----
    # Mirror extract_factors() in .ife_fit() exactly:
    #   F = sqrt(TT) * first r eigenvectors of W W' / (N*TT)   when TT <= N
    #   F = W %*% Lam / N  where Lam from N x N version        when TT > N
    if (TT <= N) {
      M     <- tcrossprod(W_fill) / (N * TT)              # TT x TT
      eig   <- eigen(M, symmetric = TRUE)
      F_new <- eig$vectors[, seq_len(r), drop = FALSE] * sqrt(TT)  # TT x r
    } else {
      M       <- crossprod(W_fill) / (N * TT)             # N x N
      eig     <- eigen(M, symmetric = TRUE)
      Lam_tmp <- eig$vectors[, seq_len(r), drop = FALSE] * sqrt(N) # N x r
      F_new   <- W_fill %*% Lam_tmp / N                            # TT x r
    }

    # ---- Lambda update: Lambda = W_fill' F / TT ----
    Lambda_new <- crossprod(W_fill, F_new) / TT            # N x r

    # ---- Convergence: check change in fitted values at observed cells ----
    # NOTE: compare fitted values, NOT F directly —
    # F is identified only up to rotation, so |F_new - F_prev| is meaningless.
    fitted_new  <- F_new      %*% t(Lambda_new)            # TT x N
    fitted_prev <- F_prev     %*% t(Lambda_prev)           # TT x N
    delta <- max(abs((fitted_new - fitted_prev)[obs_mask]))

    # Update for next iteration
    F_prev      <- F_new
    Lambda_prev <- Lambda_new

    # E-step: re-impute missing cells with updated fitted values
    W_fill[!obs_mask] <- fitted_new[!obs_mask]

    if (delta < tol_em) {
      return(list(
        F_hat       = F_new,
        Lambda_hat  = Lambda_new,
        n_iter_em   = h,
        converged_em = TRUE
      ))
    }
  }

  # Did not converge within max_iter_em
  list(
    F_hat        = F_new,
    Lambda_hat   = Lambda_new,
    n_iter_em    = max_iter_em,
    converged_em = FALSE
  )
}


# ----------------------------------------------------------------------------
# .ife_fit_unb — Outer alternating loop for unbalanced IFE
#
# Alternates between:
#   (a) EM inner loop to update F and Lambda given beta
#   (b) Unit-specific FWL projection to update beta given F
#
# @param Y_long    n_obs numeric vector of observed outcomes
# @param X_long    n_obs x p matrix of observed covariates (0 cols if p=0)
# @param unit_idx  n_obs integer vector: unit index (1-based) for each obs
# @param time_idx  n_obs integer vector: time index (1-based) for each obs
# @param N         number of unique units
# @param TT        number of unique time periods (= max T_i)
# @param r         number of factors (>= 1)
# @param tol        outer-loop convergence tolerance on max|beta_new - beta|
# @param max_iter   maximum outer-loop iterations
# @param tol_em     inner EM convergence tolerance
# @param max_iter_em maximum inner EM iterations per outer step
#
# @return list with beta, F_hat, Lambda_hat, X_tilde_long, u_long,
#              obs_mask, obs_lin, unit_rows, unit_obs_t, n_iter, converged
# ----------------------------------------------------------------------------
.ife_fit_unb <- function(Y_long, X_long, unit_idx, time_idx,
                         N, TT, r,
                         tol         = 1e-9,
                         max_iter    = 10000L,
                         tol_em      = 1e-7,
                         max_iter_em = 500L) {

  n_obs <- length(Y_long)
  p     <- ncol(X_long)

  # ---- Precompute obs_mask and linear indices into T x N matrix ----
  # Linear index (column-major): obs_lin[k] = (unit_idx[k]-1)*TT + time_idx[k]
  obs_lin  <- (unit_idx - 1L) * TT + time_idx
  obs_mask <- matrix(FALSE, TT, N)
  obs_mask[obs_lin] <- TRUE

  # ---- Precompute per-unit row indices into long vectors ----
  # unit_rows[[i]] = which(unit_idx == i)  — done once, reused every outer iter
  unit_rows  <- lapply(seq_len(N), function(i) which(unit_idx == i))
  unit_obs_t <- lapply(unit_rows,  function(idx) time_idx[idx])

  # ---- r = 0 special case: plain OLS, no EM needed ----
  if (r == 0L) {
    if (p > 0L) {
      beta <- as.vector(solve(crossprod(X_long), crossprod(X_long, Y_long)))
    } else {
      beta <- numeric(0)
    }
    F_hat        <- matrix(0, TT, 0)
    Lambda_hat   <- matrix(0, N,  0)
    X_tilde_long <- X_long
    fitted_lf    <- numeric(n_obs)   # all zero when r = 0
    u_long <- Y_long
    if (p > 0L) u_long <- u_long - as.vector(X_long %*% beta)
    return(list(
      beta         = beta,
      F_hat        = F_hat,
      Lambda_hat   = Lambda_hat,
      X_tilde_long = X_tilde_long,
      u_long       = u_long,
      obs_mask     = obs_mask,
      obs_lin      = obs_lin,
      unit_rows    = unit_rows,
      unit_obs_t   = unit_obs_t,
      n_iter       = 0L,
      converged    = TRUE
    ))
  }

  # ---- Initialisation ----
  if (p > 0L) {
    beta <- as.vector(solve(crossprod(X_long), crossprod(X_long, Y_long)))
  } else {
    beta <- numeric(0)
  }
  F_hat        <- matrix(0, TT, r)
  Lambda_hat   <- matrix(0, N,  r)
  X_tilde_long <- X_long   # will be overwritten unit-by-unit

  converged <- FALSE
  n_iter    <- 0L

  # ---- Outer loop ----
  for (iter in seq_len(max_iter)) {
    n_iter   <- iter
    beta_old <- beta

    # (a) Residual from current beta at observed cells
    if (p > 0L) {
      W_long <- Y_long - as.vector(X_long %*% beta)
    } else {
      W_long <- Y_long
    }

    # Scatter W into T x N matrix (missing cells = 0; EM will impute)
    W_mat <- matrix(0, TT, N)
    W_mat[obs_lin] <- W_long

    # (b) Inner EM loop — warm-started from current F_hat / Lambda_hat
    em <- .ife_em_inner(
      W_mat       = W_mat,
      obs_mask    = obs_mask,
      Lambda_prev = Lambda_hat,
      F_prev      = F_hat,
      N = N, TT = TT, r = r,
      tol_em = tol_em, max_iter_em = max_iter_em
    )
    F_hat      <- em$F_hat
    Lambda_hat <- em$Lambda_hat

    if (!em$converged_em) {
      warning("EM inner loop did not converge at outer iteration ", iter,
              ". Consider increasing max_iter_em or tol_em.")
    }

    # (c) Unit-specific FWL projection and beta update
    if (p > 0L) {
      A_mat <- matrix(0, p, p)
      b_vec <- numeric(p)

      for (i in seq_len(N)) {
        obs_i <- unit_rows[[i]]    # indices into long vectors for unit i
        t_i   <- unit_obs_t[[i]]   # observed time indices for unit i
        T_i   <- length(t_i)

        if (T_i == 0L) next

        F_i <- F_hat[t_i, , drop = FALSE]   # T_i x r
        X_i <- X_long[obs_i, , drop = FALSE] # T_i x p
        Y_i <- Y_long[obs_i]                  # T_i

        # FWL: X_tilde_i = (I - F_i (F_i'F_i)^{-1} F_i') X_i
        FiFi_inv  <- solve(crossprod(F_i))                        # r x r
        FtX_i     <- crossprod(F_i, X_i)                          # r x p
        FtY_i     <- crossprod(F_i, Y_i)                          # r x 1
        X_tilde_i <- X_i - F_i %*% (FiFi_inv %*% FtX_i)          # T_i x p
        Y_tilde_i <- Y_i - as.vector(F_i %*% (FiFi_inv %*% FtY_i)) # T_i

        A_mat <- A_mat + crossprod(X_tilde_i)
        b_vec <- b_vec + as.vector(crossprod(X_tilde_i, Y_tilde_i))

        X_tilde_long[obs_i, ] <- X_tilde_i   # store projected X for SE
      }

      # (d) Update beta
      beta_new <- as.vector(solve(A_mat, b_vec))

      # (e) Convergence check (outer loop)
      if (max(abs(beta_new - beta_old)) < tol) {
        beta      <- beta_new
        converged <- TRUE
        break
      }
      beta <- beta_new

    } else {
      # p = 0: nothing to iterate on; factors converge via EM alone
      X_tilde_long <- NULL
      converged    <- TRUE
      break
    }
  }

  # ---- Final residuals at observed cells ----
  # u_it = Y_it - X_it' beta - lambda_i' F_t
  fitted_lf <- rowSums(Lambda_hat[unit_idx, , drop = FALSE] *
                       F_hat[time_idx,      , drop = FALSE])  # n_obs
  u_long <- Y_long - fitted_lf
  if (p > 0L) u_long <- u_long - as.vector(X_long %*% beta)

  list(
    beta         = beta,
    F_hat        = F_hat,
    Lambda_hat   = Lambda_hat,
    X_tilde_long = X_tilde_long,
    u_long       = u_long,
    obs_mask     = obs_mask,
    obs_lin      = obs_lin,
    unit_rows    = unit_rows,
    unit_obs_t   = unit_obs_t,
    n_iter       = n_iter,
    converged    = converged
  )
}


# ----------------------------------------------------------------------------
# .ife_se_unb — Sandwich SE for unbalanced IFE
#
# Computes Var(beta_hat) = A^{-1} B A^{-1} where A and B are formed from
# projected covariates X_tilde and residuals u, summed over observed cells.
#
# @param beta          p-vector of estimated coefficients
# @param X_tilde_long  n_obs x p matrix of FWL-projected covariates
# @param u_long        n_obs residuals at observed cells
# @param unit_idx      n_obs integer unit indices
# @param N             number of units
# @param TT            number of time periods
# @param r             number of factors
# @param se_type       "standard" | "robust" | "cluster"
# @param n_obs         total observed cells
#
# @return list(vcov_mat, df)
# ----------------------------------------------------------------------------
.ife_se_unb <- function(beta, X_tilde_long, u_long, unit_idx,
                        N, TT, r, se_type, n_obs) {

  p <- length(beta)
  if (p == 0L) stop("Cannot compute SE with no covariates.")

  # ---- Degrees of freedom ----
  # Params: p (beta) + r*(N + TT - r) (interactive FE, Bai 2009 Thm 1)
  k_total <- p + r * (N + TT - r)
  df <- n_obs - k_total
  if (df <= 0L)
    stop("Degrees of freedom = ", df, " <= 0. Reduce r or use a larger panel.")

  # ---- A matrix: X_tilde' X_tilde (p x p) ----
  A     <- crossprod(X_tilde_long)   # p x p
  A_inv <- solve(A)

  # ---- Variance estimator ----
  if (se_type == "standard") {

    sigma2   <- sum(u_long^2) / df
    vcov_mat <- sigma2 * A_inv

  } else if (se_type == "robust") {

    # HC1: B = X_tilde' diag(u^2) X_tilde
    # Vectorised: scale each row of X_tilde by u_it, then crossprod
    B       <- crossprod(X_tilde_long * u_long)   # p x p
    corr    <- n_obs / (n_obs - p)                 # HC1 finite-sample correction
    vcov_mat <- A_inv %*% B %*% A_inv * corr

  } else {   # "cluster"

    # Cluster by unit i: B = sum_i psi_i psi_i', psi_i = sum_{t in T_i} u_it x_tilde_it
    B <- matrix(0, p, p)
    for (i in seq_len(N)) {
      obs_i <- which(unit_idx == i)
      if (length(obs_i) == 0L) next
      Xt_i  <- X_tilde_long[obs_i, , drop = FALSE]   # T_i x p
      u_i   <- u_long[obs_i]                           # T_i
      psi_i <- as.vector(crossprod(Xt_i, u_i))         # p
      B     <- B + tcrossprod(psi_i)                    # p x p outer product
    }
    corr     <- (N / (N - 1L)) * ((n_obs - 1L) / (n_obs - p))
    vcov_mat <- A_inv %*% B %*% A_inv * corr

  }

  list(vcov_mat = vcov_mat, df = df)
}


# ----------------------------------------------------------------------------
# ife_unbalanced — Public wrapper
#
# Fits the unbalanced panel IFE model via the Bai (2009) Appendix B EM
# algorithm. No additive fixed effects are included.
#
#' Unbalanced Panel Interactive Fixed Effects Estimator
#'
#' Fits the pure interactive fixed effects model
#' \deqn{Y_{it} = X_{it}'\beta + \lambda_i'F_t + u_{it}}
#' for unbalanced panels (units observed at different sets of time periods)
#' using the expectation-maximisation (EM) algorithm of
#' Bai (2009) Appendix B.  No additive unit or time fixed effects are
#' included in this specification.
#'
#' @param formula R formula: \code{outcome ~ covariate1 + covariate2 + ...}
#' @param data    Data frame in long format (one row per observed unit-time
#'   pair).
#' @param index   Character vector of length 2: \code{c("unit_col", "time_col")}.
#' @param r       Positive integer. Number of interactive factors (default 1).
#' @param se      SE type: \code{"standard"} (homoskedastic),
#'   \code{"robust"} (HC1), or \code{"cluster"} (cluster-robust by unit).
#'   Default \code{"standard"}.
#' @param tol     Outer-loop convergence tolerance on
#'   \eqn{\max|\hat\beta^{new} - \hat\beta^{old}|}. Default \code{1e-9}.
#' @param max_iter Maximum outer-loop iterations. Default \code{10000L}.
#' @param tol_em  Inner EM convergence tolerance on the maximum change in
#'   fitted values \eqn{\hat\lambda_i'\hat F_t} at observed cells.
#'   Default \code{1e-7}.
#' @param max_iter_em Maximum inner EM iterations per outer step.
#'   Default \code{500L}.
#'
#' @return An S3 object of class \code{"ife_unb"} with components:
#' \describe{
#'   \item{coef}{Named p-vector of estimated coefficients \eqn{\hat\beta}.}
#'   \item{vcov}{p x p variance-covariance matrix.}
#'   \item{se}{Named p-vector of standard errors.}
#'   \item{tstat}{Named p-vector of t-statistics.}
#'   \item{pval}{Named p-vector of two-sided p-values.}
#'   \item{ci}{p x 2 matrix of 95% confidence intervals.}
#'   \item{table}{Data frame coefficient table.}
#'   \item{F_hat}{TT x r estimated factor matrix (normalised F'F/TT = I_r).}
#'   \item{Lambda_hat}{N x r estimated loading matrix.}
#'   \item{residuals}{n_obs numeric vector of full-model residuals at
#'     observed cells.}
#'   \item{sigma2}{Estimated error variance (\eqn{sum(u^2)/df}).}
#'   \item{df}{Residual degrees of freedom.}
#'   \item{n_obs}{Number of observed unit-time cells.}
#'   \item{n_iter}{Outer-loop iterations to convergence.}
#'   \item{converged}{Logical.}
#'   \item{N, TT, r, se_type}{Model dimensions and options.}
#'   \item{y_name, x_names, id_col, time_col}{Variable names.}
#'   \item{unit_vals, time_vals}{Unique unit and time identifiers.}
#'   \item{unit_idx, time_idx}{Integer index vectors for \code{residuals}.}
#'   \item{call}{Matched call.}
#' }
#'
#' @references
#' Bai, J. (2009). Panel data models with interactive fixed effects.
#' \emph{Econometrica}, 77(4), 1229--1279. \doi{10.3982/ECTA6135}
#'
#' Stock, J.H. and Watson, M.W. (1998). Diffusion indexes.
#' \emph{NBER Working Paper} 6702.
#'
#' @importFrom stats pt qt
#' @export
#'
#' @examples
#' data(cigar, package = "xtife")
#' # Drop ~10 % of rows to create an unbalanced panel
#' set.seed(1)
#' cigar_unb <- cigar[sample(nrow(cigar), 1200L), ]
#' fit <- ife_unbalanced(sales ~ price, data = cigar_unb,
#'                       index = c("state", "year"), r = 2L)
#' print(fit)
# ----------------------------------------------------------------------------
ife_unbalanced <- function(formula,
                           data,
                           index,
                           r           = 1L,
                           se          = "standard",
                           tol         = 1e-9,
                           max_iter    = 10000L,
                           tol_em      = 1e-7,
                           max_iter_em = 500L) {

  cl <- match.call()

  # ================================================================
  # Input validation (V1 – V12)
  # ================================================================

  # V1 — formula
  if (!inherits(formula, "formula"))
    stop("'formula' must be an R formula object.")

  # V2 — data
  if (!is.data.frame(data))
    stop("'data' must be a data.frame.")

  # V3 — index
  if (!is.character(index) || length(index) != 2L)
    stop("'index' must be a character vector of length 2: c('unit_col', 'time_col').")

  # V4 — index columns exist
  missing_idx <- setdiff(index, names(data))
  if (length(missing_idx) > 0L)
    stop("'index' column(s) not found in 'data': ",
         paste(missing_idx, collapse = ", "))

  # V5 — se
  if (!se %in% c("standard", "robust", "cluster"))
    stop("'se' must be one of: 'standard', 'robust', 'cluster'.")

  # V6 — r
  r <- as.integer(r)
  if (r < 1L)
    stop("'r' must be a positive integer (>= 1). ",
         "For r = 0 (plain OLS), use lm() or plm().")

  # V7 — tolerances
  tol         <- as.double(tol)
  tol_em      <- as.double(tol_em)
  max_iter    <- as.integer(max_iter)
  max_iter_em <- as.integer(max_iter_em)
  if (tol    <= 0) stop("'tol' must be a positive number.")
  if (tol_em <= 0) stop("'tol_em' must be a positive number.")
  if (max_iter    < 1L) stop("'max_iter' must be >= 1.")
  if (max_iter_em < 1L) stop("'max_iter_em' must be >= 1.")

  # V8 — parse formula and check variable names
  vars    <- all.vars(formula)
  y_name  <- vars[1L]
  x_names <- vars[-1L]
  p       <- length(x_names)
  id_col   <- index[1L]
  time_col <- index[2L]

  all_needed <- c(y_name, x_names, id_col, time_col)
  missing_v  <- setdiff(all_needed, names(data))
  if (length(missing_v) > 0L)
    stop("Variable(s) not found in 'data': ",
         paste(missing_v, collapse = ", "))

  # V9 — no duplicate (unit, time) pairs
  dup_key <- paste(data[[id_col]], data[[time_col]], sep = "___IFE___")
  if (anyDuplicated(dup_key))
    stop("Duplicate (unit, time) pairs found. ",
         "Each (i, t) must appear at most once in 'data'.")

  # V10 — no NA in outcome or covariates
  for (v in c(y_name, x_names)) {
    if (anyNA(data[[v]]))
      stop("Missing values (NA) found in variable '", v, "'. ",
           "Structural missing observations should be represented by ",
           "absent rows, not NA values.")
  }

  # ---- Sort by unit then time (required for index alignment) ----
  data <- data[order(data[[id_col]], data[[time_col]]), ]

  unit_vals <- unique(data[[id_col]])
  time_vals <- unique(data[[time_col]])
  N  <- length(unit_vals)
  TT <- length(time_vals)

  unit_idx <- match(data[[id_col]],   unit_vals)   # integer 1..N
  time_idx <- match(data[[time_col]], time_vals)   # integer 1..TT
  n_obs    <- nrow(data)

  # V11 — each unit needs at least r+1 observations (FWL requires T_i > r)
  obs_per_unit <- tabulate(unit_idx, nbins = N)
  if (any(obs_per_unit < r + 1L))
    stop("Some units have fewer than r + 1 = ", r + 1L,
         " observations. The FWL projection requires T_i > r for all units. ",
         "Reduce r or remove units with too few observations.")

  # V12 — r vs dimensions
  if (r > min(N, TT))
    stop("r = ", r, " exceeds min(N, TT) = min(", N, ", ", TT, ") = ",
         min(N, TT), ". Reduce r.")

  # ================================================================
  # Build long vectors
  # ================================================================
  Y_long <- as.double(data[[y_name]])
  X_long <- if (p > 0L) {
    as.matrix(data[, x_names, drop = FALSE])
  } else {
    matrix(0, n_obs, 0L)
  }
  storage.mode(X_long) <- "double"

  # ================================================================
  # Grand-mean centering (matches ife(force = "none") behaviour)
  #
  # Subtracting the grand mean from Y and X:
  #   (a) gives the same OLS initialisation as grand-mean-demeaned OLS
  #       (avoids the "regression-through-origin" trap when Y has a large
  #        mean, e.g. sales ≈ 100 for the cigar panel)
  #   (b) ensures the algorithm converges to the same solution as
  #       ife(force = "none") on balanced panels
  #   (c) is valid for unbalanced panels (grand mean over observed cells)
  #
  # Note: beta is invariant to grand-mean centering when r >= 1 factors
  # can absorb the constant direction (F_t ≈ 1 ⟹ lambda_i absorbs offset).
  # ================================================================
  mu_Y <- mean(Y_long)
  Y_long_c <- Y_long - mu_Y

  if (p > 0L) {
    mu_X     <- colMeans(X_long)                                    # p-vector
    X_long_c <- X_long - matrix(mu_X, nrow = n_obs, ncol = p,
                                 byrow = TRUE)
  } else {
    mu_X     <- numeric(0L)
    X_long_c <- X_long                                              # 0-col matrix
  }

  # ================================================================
  # Core estimation (on centred data)
  # ================================================================
  fit <- .ife_fit_unb(
    Y_long      = Y_long_c,
    X_long      = X_long_c,
    unit_idx    = unit_idx,
    time_idx    = time_idx,
    N = N, TT = TT, r = r,
    tol         = tol,
    max_iter    = max_iter,
    tol_em      = tol_em,
    max_iter_em = max_iter_em
  )

  if (!fit$converged)
    warning("Outer loop did not converge after ", max_iter,
            " iterations. Increase max_iter or relax tol.")

  # ================================================================
  # Standard errors
  # ================================================================
  coef_vec <- fit$beta
  names(coef_vec) <- x_names

  vcov_mat <- matrix(NA_real_, p, p,
                     dimnames = list(x_names, x_names))
  se_vec   <- rep(NA_real_, p);  names(se_vec) <- x_names
  tstat    <- rep(NA_real_, p);  names(tstat)  <- x_names
  pval     <- rep(NA_real_, p);  names(pval)   <- x_names
  ci       <- matrix(NA_real_, p, 2L,
                     dimnames = list(x_names, c("2.5 %", "97.5 %")))
  df_resid <- n_obs - p - r * (N + TT - r)

  if (p > 0L) {
    se_list  <- .ife_se_unb(
      beta         = fit$beta,
      X_tilde_long = fit$X_tilde_long,
      u_long       = fit$u_long,
      unit_idx     = unit_idx,
      N = N, TT = TT, r = r,
      se_type = se,
      n_obs   = n_obs
    )
    vcov_mat  <- se_list$vcov_mat
    df_resid  <- se_list$df
    se_vec    <- sqrt(pmax(diag(vcov_mat), 0))
    names(se_vec) <- x_names
    tstat     <- coef_vec / se_vec
    pval      <- 2 * pt(-abs(tstat), df = df_resid)
    t_crit    <- qt(0.975, df = df_resid)
    ci[, 1L]  <- coef_vec - t_crit * se_vec
    ci[, 2L]  <- coef_vec + t_crit * se_vec
  }

  sigma2 <- if (df_resid > 0L) sum(fit$u_long^2) / df_resid else NA_real_

  # ================================================================
  # Coefficient table (mirroring print.ife format)
  # ================================================================
  coef_table <- if (p > 0L) {
    data.frame(
      Estimate  = coef_vec,
      Std.Error = se_vec,
      t.value   = tstat,
      Pr.t      = pval,
      CI.lower  = ci[, 1L],
      CI.upper  = ci[, 2L],
      row.names = x_names,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }

  # ================================================================
  # Return
  # ================================================================
  structure(
    list(
      coef        = coef_vec,
      vcov        = vcov_mat,
      se          = se_vec,
      tstat       = tstat,
      pval        = pval,
      ci          = ci,
      table       = coef_table,
      F_hat       = fit$F_hat,
      Lambda_hat  = fit$Lambda_hat,
      residuals   = fit$u_long,
      sigma2      = sigma2,
      df          = df_resid,
      n_obs       = n_obs,
      n_iter      = fit$n_iter,
      converged   = fit$converged,
      N           = N,
      TT          = TT,
      r           = r,
      se_type     = se,
      y_name      = y_name,
      x_names     = x_names,
      id_col      = id_col,
      time_col    = time_col,
      unit_vals   = unit_vals,
      time_vals   = time_vals,
      unit_idx    = unit_idx,
      time_idx    = time_idx,
      call        = cl
    ),
    class = "ife_unb"
  )
}


# ----------------------------------------------------------------------------
# print.ife_unb — S3 print method
#
#' @export
# ----------------------------------------------------------------------------
print.ife_unb <- function(x, digits = 4L, ...) {

  cat("\n")
  cat("Unbalanced Panel IFE  (Bai 2009 Appendix B EM Algorithm)\n")
  cat(strrep("-", 58L), "\n")
  cat(sprintf("Panel    : N = %d units,  TT = %d periods (max)\n",
              x$N, x$TT))
  cat(sprintf("Observed : n_obs = %d cells  (%.1f%% of N x TT)\n",
              x$n_obs, 100 * x$n_obs / (x$N * x$TT)))
  cat(sprintf("Factors  : r = %d\n", x$r))
  se_label <- switch(x$se_type,
    "standard" = "standard (homoskedastic)",
    "robust"   = "robust (HC1)",
    "cluster"  = paste0("cluster-robust by unit (", x$id_col, ")")
  )
  cat(sprintf("SE type  : %s\n", se_label))
  cat(sprintf("Outcome  : %s\n", x$y_name))
  cat(strrep("-", 58L), "\n")

  if (nrow(x$table) > 0L) {
    tbl   <- x$table
    stars <- ifelse(tbl$Pr.t < 0.01, "***",
              ifelse(tbl$Pr.t < 0.05, "**",
               ifelse(tbl$Pr.t < 0.10, "*", "")))
    out <- data.frame(
      Estimate   = formatC(tbl$Estimate,  digits = digits, format = "f"),
      Std.Error  = formatC(tbl$Std.Error, digits = digits, format = "f"),
      t.value    = formatC(tbl$t.value,   digits = digits, format = "f"),
      "Pr(>|t|)" = formatC(tbl$Pr.t,      digits = digits, format = "f"),
      "95% CI"   = paste0("[",
                     formatC(tbl$CI.lower, digits = digits, format = "f"),
                     ", ",
                     formatC(tbl$CI.upper, digits = digits, format = "f"),
                     "]"),
      " "        = stars,
      row.names  = rownames(tbl),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    print(out, quote = FALSE)
    cat("---\n")
    cat("Signif. codes: *** < 0.01  ** < 0.05  * < 0.10\n")
  } else {
    cat("(No covariates specified)\n")
  }

  cat(strrep("-", 58L), "\n")
  cat(sprintf("sigma^2 = %.6g  |  df = %d\n", x$sigma2, x$df))
  cat(sprintf("Converged: %s  |  Outer iterations: %d\n",
              if (x$converged) "YES" else "NO (increase max_iter)",
              x$n_iter))
  cat("\n")
  invisible(x)
}
