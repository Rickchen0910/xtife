# ==============================================================================
# collate_mc.R
# Reads all per-cell .rds files from results/ and assembles JSS Tables 1–5
# as both plain R data frames and LaTeX (via xtable).
#
# Run locally after downloading results/ from HPC:
#   rsync -av bc25911@hpc.essex.ac.uk:/home/bc25911/xtife_mc/results/ ./results/
#   Rscript collate_mc.R
#
# Output (written to tables/ next to this script):
#   tab1_static_bal.tex   tab2_static_unb_f80.tex   tab2_static_unb_f60.tex
#   tab3_dyn_bal.tex      tab4_dyn_unb_f80.tex      tab4_dyn_unb_f60.tex
#   tab5_factor_sel.tex
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("xtable", quietly = TRUE))
    install.packages("xtable", repos = "https://cloud.r-project.org")
  library(xtable)
})

# ==============================================================================
# 0 — Paths
# ==============================================================================

## Locate this script's directory whether run interactively or via Rscript
if (interactive()) {
  script_dir <- getwd()
} else {
  argv       <- commandArgs(trailingOnly = FALSE)
  file_arg   <- argv[grep("^--file=", argv)][1L]
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg),
                                      mustWork = FALSE))
}
results_dir <- file.path(script_dir, "results")
tables_dir  <- file.path(script_dir, "tables")
dir.create(tables_dir, showWarnings = FALSE)

rds_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0L)
  stop("No .rds files found in: ", results_dir, "\n  Run jobs first.")

cat(sprintf("Found %d result files in %s\n", length(rds_files), results_dir))

# ==============================================================================
# 1 — Parse filenames and load all results
# ==============================================================================

parse_fname <- function(f) {
  # Pattern: dgp{d}_N{N}_T{TT}_f{fill}_{B}reps.rds
  b <- sub("\\.rds$", "", basename(f))
  p <- strsplit(b, "_")[[1L]]
  list(
    dgp_id   = as.integer(sub("dgp", "", p[1L])),
    N        = as.integer(sub("N",   "", p[2L])),
    TT       = as.integer(sub("T",   "", p[3L])),
    fill_int = as.integer(sub("f",   "", p[4L])),
    B        = as.integer(sub("reps","", p[5L]))
  )
}

all_rows <- lapply(rds_files, function(f) {
  tryCatch({
    obj    <- readRDS(f)
    params <- parse_fname(f)
    row    <- data.frame(
      dgp_id   = params$dgp_id,
      N        = params$N,
      TT       = params$TT,
      fill     = params$fill_int / 100,
      fill_int = params$fill_int,
      B        = params$B,
      stringsAsFactors = FALSE
    )
    cbind(row, obj$summ)
  }, error = function(e) {
    warning("Could not read: ", basename(f), " — ", conditionMessage(e))
    NULL
  })
})

df_all <- do.call(rbind, Filter(Negate(is.null), all_rows))
df_all <- df_all[order(df_all$dgp_id, df_all$fill_int, df_all$N, df_all$TT), ]
cat(sprintf("Loaded %d cells total.\n", nrow(df_all)))

# ==============================================================================
# 2 — Formatting helpers
# ==============================================================================

fmt  <- function(x, d = 3L) formatC(x, digits = d, format = "f")
fmtp <- function(x, d = 3L) {                   # for coverage: 0.950 → ".950"
  s <- formatC(x, digits = d, format = "f")
  sub("^0", "", s)
}

## Build one rows × stats data frame for a given subset
make_tab_df <- function(sub) {
  data.frame(
    N          = sub$N,
    T          = sub$TT,
    OLS.Bias   = fmt(sub$OLS.Bias),
    IFE.Bias   = fmt(sub$IFE.Bias),
    IFE.SD     = fmt(sub$IFE.SD),
    IFE.RMSE   = fmt(sub$IFE.RMSE),
    BC.Bias    = fmt(sub$BC.Bias),
    BC.RMSE    = fmt(sub$BC.RMSE),
    SE.SD      = fmt(sub$IFE.SE_SD),
    Cov.Std    = fmtp(sub$IFE.Cov.Std),
    Cov.Rob    = fmtp(sub$IFE.Cov.Rob),
    Cov.Cl     = fmtp(sub$IFE.Cov.Cl),
    BC.Cov     = fmtp(sub$BC.Cov.Std),
    Size       = fmtp(sub$IFE.Size),
    r.pct      = fmt(sub$r.correct * 100, 1L),   # % correct r selection
    Conv       = fmt(sub$conv.rate, 3L),
    stringsAsFactors = FALSE
  )
}

## Column labels for xtable
COL_NAMES <- c("$N$", "$T$",
               "OLS", "IFE", "SD", "RMSE",
               "BC", "BC.RMSE",
               "SE/SD",
               "Std", "Rob", "Cl", "BC.Std",
               "Size",
               "$P(\\hat r=r)$", "Conv.")

## Write one .tex table file
write_tex <- function(tab_df, file, caption, label,
                      add_header = TRUE) {
  xt <- xtable(tab_df,
               caption = caption,
               label   = label,
               digits  = 0L)             # all pre-formatted as strings
  names(xt) <- COL_NAMES

  ## Multicolumn header for bias / coverage groups
  addtorow <- list()
  if (add_header) {
    addtorow$pos <- list(-1L)
    addtorow$command <- paste0(
      "\\hline\n",
      " & & \\multicolumn{6}{c}{Bias \\& RMSE} & ",
      "\\multicolumn{5}{c}{95\\% Coverage} & & & \\\\\n",
      "\\cline{3-8}\\cline{9-13}\n"
    )
  }

  print(xt,
        file             = file,
        include.rownames = FALSE,
        include.colnames = TRUE,
        add.to.row       = if (add_header) addtorow else NULL,
        sanitize.text.function = identity,   # keep LaTeX commands in colnames
        hline.after      = c(-1L, 0L, nrow(tab_df)),
        floating         = TRUE,
        caption.placement = "top")

  cat(sprintf("  Written: %s\n", file))
}

# ==============================================================================
# 3 — Build and write each table
# ==============================================================================

cat("\nBuilding tables ...\n")

## --- Table 1: Static balanced (DGP 1) ---------------------------------------
sub1 <- df_all[df_all$dgp_id == 1L, ]
if (nrow(sub1) > 0L) {
  tab1 <- make_tab_df(sub1)
  write_tex(tab1,
    file    = file.path(tables_dir, "tab1_static_bal.tex"),
    caption = paste("Monte Carlo results: static balanced panel (DGP 1).",
                    "\\emph{N}: units; \\emph{T}: periods;",
                    "OLS/IFE/BC: bias of each estimator;",
                    "SE/SD: ratio of mean SE to Monte Carlo SD;",
                    "Cov: empirical 95\\% coverage; Size: rejection rate under $H_0$."),
    label   = "tab:mc_static_bal")
} else {
  cat("  [SKIP] No DGP 1 results found.\n")
}

## --- Tables 2a/2b: Static unbalanced (DGP 2) --------------------------------
for (fi in sort(unique(df_all$fill_int[df_all$dgp_id == 2L]))) {
  sub2 <- df_all[df_all$dgp_id == 2L & df_all$fill_int == fi, ]
  if (nrow(sub2) == 0L) next
  tab2 <- make_tab_df(sub2)
  write_tex(tab2,
    file    = file.path(tables_dir, sprintf("tab2_static_unb_f%d.tex", fi)),
    caption = sprintf(
      "Monte Carlo results: static unbalanced panel (DGP 2, fill = %d\\%%).",
      fi),
    label   = sprintf("tab:mc_static_unb_f%d", fi))
}

## --- Table 3: Dynamic balanced (DGP 3) --------------------------------------
sub3 <- df_all[df_all$dgp_id == 3L, ]
if (nrow(sub3) > 0L) {
  ## For dynamic balanced: replace Cov.Cl column with Cov.HAC label
  tab3 <- make_tab_df(sub3)
  names(tab3)[names(tab3) == "Cov.Cl"] <- "Cov.Cl*"  # cluster = HAC proxy here
  write_tex(tab3,
    file    = file.path(tables_dir, "tab3_dyn_bal.tex"),
    caption = paste("Monte Carlo results: dynamic balanced panel (DGP 3).",
                    "Regressors: $x_{it} = y_{i,t-1}$; bias correction via Moon \\&",
                    "Weidner (2017). $T \\geq 50$ only (burn-in requirement)."),
    label   = "tab:mc_dyn_bal")
} else {
  cat("  [SKIP] No DGP 3 results found.\n")
}

## --- Tables 4a/4b: Dynamic unbalanced (DGP 4) --------------------------------
for (fi in sort(unique(df_all$fill_int[df_all$dgp_id == 4L]))) {
  sub4 <- df_all[df_all$dgp_id == 4L & df_all$fill_int == fi, ]
  if (nrow(sub4) == 0L) next
  tab4 <- make_tab_df(sub4)
  ## Replace Cov.Cl with Cov.HAC for DGP 4 (HAC SE is the correct one here)
  tab4$Cov.Cl <- fmtp(sub4$IFE.Cov.HAC)
  names(tab4)[names(tab4) == "Cov.Cl"] <- "Cov.HAC"
  write_tex(tab4,
    file    = file.path(tables_dir, sprintf("tab4_dyn_unb_f%d.tex", fi)),
    caption = sprintf(
      "Monte Carlo results: dynamic unbalanced panel (DGP 4, fill = %d\\%%).",
      fi),
    label   = sprintf("tab:mc_dyn_unb_f%d", fi))
}

## --- Table 5: Factor selection accuracy (all DGPs) ---------------------------
if (nrow(df_all) > 0L) {
  tab5 <- data.frame(
    DGP      = df_all$dgp_id,
    Fill     = paste0(df_all$fill_int, "\\%"),
    N        = df_all$N,
    T        = df_all$TT,
    r.mean   = fmt(df_all$r.mean, 2L),
    r.pct    = fmt(df_all$r.correct * 100, 1L),
    r.under  = fmt(df_all$r.under  * 100, 1L),
    r.over   = fmt(df_all$r.over   * 100, 1L),
    stringsAsFactors = FALSE
  )
  xt5 <- xtable(tab5,
    caption = "Factor selection accuracy across all DGPs ($r_{\\rm true} = 2$, $r_{\\max} = 5$).",
    label   = "tab:mc_factor_sel",
    digits  = 0L)
  names(xt5) <- c("DGP", "Fill", "$N$", "$T$",
                   "$\\bar r$", "$P(\\hat r=2)$\\%",
                   "$P(\\hat r<2)$\\%", "$P(\\hat r>2)$\\%")
  print(xt5,
    file             = file.path(tables_dir, "tab5_factor_sel.tex"),
    include.rownames = FALSE,
    sanitize.text.function = identity,
    hline.after      = c(-1L, 0L, nrow(tab5)),
    floating         = TRUE,
    caption.placement = "top")
  cat(sprintf("  Written: %s\n", file.path(tables_dir, "tab5_factor_sel.tex")))
}

# ==============================================================================
# 4 — Console summary: check key statistics against JSS criteria
# ==============================================================================

cat("\n=== Key statistic check (JSS reviewer targets) ===\n")
cat(sprintf("%-30s  %8s  %8s\n", "Metric", "Observed", "Target"))
cat(strrep("-", 50L), "\n")

check <- function(label, val, target_lo, target_hi, val_str = NULL) {
  vs <- if (!is.null(val_str)) val_str else sprintf("%.4f", val)
  ok <- if (is.null(val)) "NO DATA" else if (val >= target_lo & val <= target_hi) "OK" else "CHECK"
  cat(sprintf("%-30s  %8s  [%.3f,%.3f]  %s\n", label, vs, target_lo, target_hi, ok))
}

## Check DGP1 (N=100, T=100 cell if available)
r100 <- df_all[df_all$dgp_id == 1L & df_all$N == 100 & df_all$TT == 100, ]
if (nrow(r100) > 0L) {
  check("DGP1 IFE Cov.Std  (N=T=100)", r100$IFE.Cov.Std[1], 0.93, 0.97)
  check("DGP1 SE/SD ratio  (N=T=100)", r100$IFE.SE_SD[1],   0.90, 1.10)
  check("DGP1 Size          (N=T=100)", r100$IFE.Size[1],    0.03, 0.07)
  check("DGP1 r.correct %  (N=T=100)", r100$r.correct[1],   0.80, 1.00)
}

## Check DGP2 f=80 (N=100, T=100)
r200 <- df_all[df_all$dgp_id == 2L & df_all$fill_int == 80 &
               df_all$N == 100 & df_all$TT == 100, ]
if (nrow(r200) > 0L) {
  check("DGP2 IFE Cov.Std  (N=T=100,f80)", r200$IFE.Cov.Std[1], 0.93, 0.97)
}

cat(strrep("-", 50L), "\n")
cat("Done. LaTeX tables written to: ", tables_dir, "\n")

# ==============================================================================
# 5 — Optional: save combined data frame for further analysis
# ==============================================================================

saveRDS(df_all,
        file = file.path(script_dir, "mc_results_combined.rds"))
cat("Combined data frame saved: mc_results_combined.rds\n")
