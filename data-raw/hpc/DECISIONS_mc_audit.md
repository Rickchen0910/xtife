# MC Audit — Decisions & Bias-Correction Verification Log

Audit of the JSS Monte Carlo study (`mc_xtife_hpc.R`) and the bias-correction code
it exercises. See plan `~/.claude/plans/polished-greeting-barto.md`.

## Part A — Bias-correction code verified against the source papers

### A1. SWW2025 unbalanced BC — `.ife_bias_unb()` in `R/ife_unbalanced.R` — BUG FOUND & FIXED

**Symptom.** Under i.i.d. homoskedastic errors the raw unbalanced IFE is unbiased
(N=T=60, 200 reps: bias +0.0008, MCSE 0.0014), but `bias_corr=TRUE` introduced a
large spurious **negative** bias (−0.076), getting worse as the panel grew more
unbalanced. The same overcorrection appeared in the dynamic-weak case (DGP 4).

**Root cause.** SWW2025 Theorem 4.2 (paper p.19, `Unbalanced IFE.pdf`) defines
the curvature inverses with a **leading minus sign**:
  `[L̄_ff']_t = (−Σ_i d_it λ_i λ_i')⁻¹`,  `[L̄_λλ']_i = (−Σ_t d_it f_t f_t')⁻¹`.
The code stored these as the *positive*-definite inverses (`solve(Σ d λλ')`), which
is correct for the δ/ω least squares (eq. 4.3) and for `Ξ̂/Δ̂` (they use the inverse
**twice**, so the two minus signs cancel). But `b3, b4, b2` and the HAC `b4` use the
inverse **once**, so they were off by a factor of −1. For balanced homoskedastic data
`b3 = b5` analytically (verified); with the wrong sign the code computed
`b̂ = b3+b4+b5+b6 ≈ 2·b5` instead of the correct `b̂ = −b3−b4+b5+b6 ≈ 0`.

**Fix.** In `.ife_bias_unb()` apply the paper's minus sign to the single-use terms:
`b3`, `b4` (i.i.d.), `b2` (weak exog), and `b4_hac`. Leave `b5`, `b6`, `b6_hac`
unchanged (double use → sign cancels). `b_hat = b2+b3+b4+b5+b6` and the δ/ω/Ξ/Δ
machinery are untouched, so test S8 (`b_hat == sum`) still holds.

**Validation (post-fix, fixed package code).**
- Static unbalanced, iid, true β=1, N=T=50, 60 reps: bc ≈ raw ≈ 1.00 at fill ∈
  {1.0, 0.8, 0.6} (was 0.927 / 0.911 / 0.891). Spurious overcorrection gone.
- Dynamic unbalanced (SWW-weak), iid, true β=0.5, N=40: BC now removes the
  Nickell-type bias — T=30 raw 0.464 → bc 0.501; T=50 raw 0.475 → bc 0.497.
- All existing tests pass (`test-ife.R`, `test-ife-unbalanced.R`,
  `test-ife-unbalanced-inference.R`) with `NOT_CRAN=true`.

### A2. Moon–Weidner dynamic BC — `.bias_correct_mw()` in `R/ife.R` — BUG FIXED
#### Bias: was correct.  Variance: was bugged → NOW FIXED.

**FIX APPLIED & VALIDATED.** Replaced the Bai (2009) `B_hat`/`C_hat` (B2/B3) reused
inside `.bias_correct_mw` with the Moon–Weidner Definition 1 estimators built from
the factor-PROJECTED regressor (M_λ X_k for B2, M_f X_k for B3); kept B1; corrected
estimator = β + W⁻¹(B1/T + B2/N + B3/T). The static `.bias_correct` is untouched.
Validated against the MW Table-1 benchmark and in the MC pipeline:
| cell | BC/FLS std (before→after) | BC.cov (before→after) |
|---|---|---|
| r=1 N=100 T=40 | 1.48 → 0.92 (MW ref 0.94) | — |
| r=2 N=100 T=40 | 1.88 → 0.85 | — |
| DGP3 N=100 T=100 (pipeline) | — | 0.535 → **0.92** |
| DGP3 N=60 T=50 (pipeline) | — | 0.57 → **0.91** |
Bias stays corrected (e.g. N=T=100: FLS −0.010 → BC −0.003). All package tests pass;
`devtools::check()` = 0/0/0.  Original diagnosis retained below for the record.

#### (Original diagnosis)
- **Scaling/sign correct.** Verified against `Dynamic_IFE.pdf` Definition 1 and
  Theorem 4.3: B1_k = (1/N) Σ_i Σ_{t<s≤t+M} Γ((s−t)/M)[P̂f]_ts ê_it x_{k,is};
  the corrected estimator subtracts (1/T)Ŵ⁻¹B̂1 (+ Bai B2/N + B3/T). The code's
  `sum(Pf*cross_trunc)/(NT)` then `D0_inv %*% B1_vec` reproduces this exactly. The
  **bias is removed correctly** at all r (matches MW Table 1 bias numbers at r=1).
- **Variance is too large** (open bug). MW Table 1 (N=100, AR(1), ρ0=0.3) shows
  BC-FLS std ≤ FLS std (ratio 0.66→0.97): their correction *reduces/preserves*
  variance. Our implementation **inflates** it:
  | cell | BC/FLS std | MW ref |
  |---|---|---|
  | r=1, T=20 | 1.13 | ~0.85 |
  | r=1, T=40 | 1.48 | 0.94 |
  | r=2, T=40 | 1.88 | — |
  | r=3, T=40 | 1.83 | — |
  Inflation is ~2× for r=2 and does NOT vanish with N,T (stable to N=T=200; not a
  κ=√(N/T) effect — tested N=40/T=200, N=200/T=40 both ~2×). The **raw** estimator
  matches MW (FLS std ≈ MW); only the BC step diverges. `M1` (1 vs 4) makes no
  difference, so the hardcoded `M1=1` is not the cause.
- **Consequence for the MC**: DGP 3 (balanced dynamic) BC reduces bias
  (−0.010→−0.004 at N=T=100) but BC coverage COLLAPSES (≈0.54), because the CI
  uses the (correct, raw-level) SE while the BC point estimate's SD is ~2× larger.
- **Root cause — CAUGHT (MW Table 1 benchmark, term-by-term decomposition).**
  Reproduced MW Table 1 (r=1, N=100, AR(1), ρ0=0.3) and decomposed the correction
  into its B1 (dynamic), B2 (cross-section) and B3 (time) parts:
  | term | mean | SD | cor(raw,·) |
  |---|---|---|---|
  | corr_B1 (dynamic)      | +0.0090 | 0.0039 | −0.49 (cancels noise) |
  | corr_B2 (= −B2_hat/N)  | ~0      | 0.0001 | — |
  | corr_B3 (= −B3_hat/T)  | ~0      | 0.0205 | +0.03 (pure noise) |
  - **B1 alone** reproduces MW: raw+corr_B1 has SD ratio **0.91** (MW: 0.94). ✓
  - **B3 (= Bai `C_hat`) is the culprit**: `B3_hat` has SD = 0.82 (vs B2_hat 0.011),
    i.e. the *time*-heteroskedasticity term is wildly noisy.
  - **Mechanism**: `C_hat` is weighted by `omega_t = (1/N)Σ_i ê_it²` (per-time
    error variance). Its sampling noise is COMMON to all units at time t, so it
    does NOT average out when the formula sums over i — and it is amplified by the
    large lagged-dependent regressor x_it = y_{i,t-1}. `B_hat`'s `sigma2_i`
    (per-unit variance) is unit-specific and averages away, so B2 is fine. Under
    homoskedasticity B3=0 in truth, so this is pure added variance.
  - **Fix direction**: in the dynamic path (`.bias_correct_mw`), replace the Bai
    `C_hat`/`B_hat` (B2/B3) with the lower-variance Moon–Weidner Definition 1
    estimators B̂2,k = (1/T)Σ ê_it²[M_λ X_k f(f'f)⁻¹(λ'λ)⁻¹λ']_ii and
    B̂3,k = (1/N)Σ ê_it²[M_f X_k λ(λ'λ)⁻¹(f'f)⁻¹f']_tt, which use the
    factor-PROJECTED regressor (removing the large factor component of x that
    amplifies the omega noise). Keep B1 (already correct). Leave the STATIC
    `.bias_correct` untouched (DGP 1 BC.cov ≈ 0.94 is fine — static x is not large).
    Verify the fix against the MW Table 1 (r=1) benchmark: BC/FLS std should drop
    from 1.50 to ≈0.94, and bias must stay corrected.  Math change → verify before
    shipping (CRAN package code).
- **Contrast**: the SWW *unbalanced* dynamic BC (DGP 4, exog="weak") does NOT have
  this problem — BC.cov ≈ 0.88 and improves with N,T. The variance issue is
  specific to the balanced Moon–Weidner path in `ife.R`.

### A3. Bai (2009) static balanced BC — `.bias_correct()` in `R/ife.R` — OK (no change)
Under i.i.d. errors the B/N and C/T terms are ≈0 (correct — they require
heteroskedasticity). No sign bug. To be re-validated against a heteroskedastic DGP.

## Part B — DGP redesign (DONE)

### Outcome
Validated empirically (30-rep mini-MC via `run_one_rep` + `summarise_reps`) that no
MCAR DGP simultaneously delivers (a) a large *correctable static* IFE bias, (b)
reliable factor selection, and (c) good coverage at moderate N,T:
- The static incidental-parameter bias is genuinely O(1/N+1/T) and **small**; the
  Bai/SWW estimators are nearly unbiased under homoskedasticity, so BC is a small
  refinement there (correct, validated not to over-correct — Part A).
- **Heteroskedastic errors** (needed to enlarge the static bias) make IC(BIC)/SVT
  **over-select** the number of factors (r.correct collapses).
- **SWW loadings λ~N(1,1)** degrade coverage (0.6–0.7) in this pipeline.
- **Factor-correlated (Pattern-2) missingness** produces a large raw bias but BC
  does not cleanly remove it at small N,T (needs slow large-N,T validation) and the
  panels converge slowly.

### Final design implemented in `mc_xtife_hpc.R`
Clean, well-behaved, homoskedastic design (well-conditioned factor structure):
- λ_i, F_t ~ N(0, I_r); u_it ~ N(0,1) i.i.d.
- **Static (DGP 1/2)**: x_it = (λ_i + μ_i)'F_t + ε_it, μ_i ~ N(0,I_r) independent.
  The extra loading μ keeps x correlated with y's factor term (large OLS bias →
  motivates IFE) WITHOUT x being a clean factor combination (which made IC select
  r=1). Validated: OLS bias ≈ +0.4, IFE unbiased, Cov.rob ≈ 0.93, r.correct ≈ 0.9.
- **Dynamic (DGP 3/4)**: x_it = y_{i,t-1}, AR(1) factors, **θ_ma = 0** (serially
  uncorrelated errors → predetermined regressor; the previous MA(1) was the fatal
  endogeneity bug). Validated: raw IFE bias ≈ −0.02, **BC removes it** (DGP 4:
  −0.021 → +0.002), r.correct ≈ 1.0.
- **beta_true**: 1.0 (static), 0.5 (dynamic, stationary).
- Missingness: Pattern-1 MCAR at `fill`.

### Deferred (documented, not implemented)
A dramatic *static* BC showcase would need SWW Pattern-2 missingness and/or strong
heteroskedasticity plus the matching HAC bias terms, with the factor-selection and
coverage interactions resolved. This is a separate simulation-design effort.

### Validation (30-rep mini-MC, N=60)
| DGP | OLS | IFE.bias | BC.bias | Cov.std | r.correct |
|-----|-----|----------|---------|---------|-----------|
| 1 static bal | +0.39 | +0.008 | +0.009 | 0.93 | 0.87 |
| 2 static unb | +0.39 | +0.002 | +0.002 | 0.87 | 0.93 |
| 3 dyn bal    | +0.22 | −0.020 | −0.008 | 0.77\* | 1.00 |
| 4 dyn unb    | +0.23 | −0.021 | +0.002 | 0.77\* | 0.97 |
\* dynamic raw coverage is low because the raw estimator is biased; it improves
after BC and at larger N,T. Coverage figures are noisy at 30 reps / N=60.

### SWW2025 §6.1 reference (for any future Pattern-2 / heteroskedastic extension)
- **Static DGP 1**: λ~N(1,1); x1 = 1 + Σ(λ+μ)(f_t+f_{t-1}) + N(0,1), μ~N(1,1);
  errors v = (e_t+e_{t-1})/√2 with e~N(0,σ²), σ~U(0.5,1.5) (heteroskedastic MA(1),
  strictly exogenous → HAC/cluster SE).
- **Dynamic DGP 2**: x1 = y_{i,t-1}; AR(1) factors (ρ_f=0.5, σ_f=0.5); errors
  v~N(0,σ²), σ~U(0.5,1.5) — heteroskedastic but **serially uncorrelated** (no MA),
  because a lagged-y regressor + serial errors would be endogenous.
- Missing patterns: Pattern 1 = MCAR (E d_it ~ U(0.5,0.9)); **Pattern 2 = selection
  on factors** (E d_it = Φ(λ_i'f_t)) — this is where the unbalanced bias is large and
  BC is essential. MCAR yields negligible bias, so the MC must use Pattern 2 (or a
  factor-correlated missingness) to demonstrate the correction.

Decision: redesign DGPs to (i) use heteroskedastic errors everywhere; (ii) keep
serial correlation only in the static cases (strictly exogenous regressor);
(iii) set `theta_ma=0` for the dynamic cases; (iv) use factor-correlated missingness
(Pattern 2) for the unbalanced cases so the IFE bias is non-trivial.

## Part A3-bis & Part D — static BC fix + complex (heteroskedastic) MC

### A3-bis. `.bias_correct()` (static Bai) — same C-term variance bug, FIXED
The time term (Bai C_hat) used the raw regressor and amplified the common
per-time variance noise omega_t — invisible under homoskedasticity (C=0) but it
OVER-corrected under heteroskedasticity (e.g. BC.bias +0.002 -> -0.010). Replaced
B/C with the Moon-Weidner Definition-1 projected forms (M_lambda X_k, M_f X_k),
identical machinery to `.bias_correct_mw`. Under homoskedasticity B=C=0 exactly
(lambda'M_lambda=0, f'M_f=0) so the simple-MC static results are unchanged; under
heteroskedasticity it no longer over-corrects (BC.bias -> +0.001). Static draw
order in the DGP generator was preserved so the homoskedastic data is byte-identical.
devtools::check() = 0/0/0; all tests pass.

### Part D — second (complex) Monte Carlo: factor-driven heteroskedasticity
Selected via env var `MC_ERR=hetero` (default `homo` = simple MC, unchanged).
Error variance tied to the common factor term: sigma_it = 0.4 + 0.8|lambda_i'F_t|/sd,
rescaled to unit mean variance. Validated (B=150, N=T=80):
| DGP | IFE.bias | BC.bias | Cov.std | Cov.rob | Cov.cl | BC | SE/SD | r.ok |
|-----|----------|---------|---------|---------|--------|----|-------|------|
| 1 static bal | +.001 | +.001 | .95 | .95 | .95 | .97 | .92 | 1.00 |
| 2 static unb | +.003 | +.003 | .92 | .91 | .91 | .91 | .91 | 1.00 |
| 3 dyn bal | -.018 | -.005 | .58 | .71 | .73 | .84 | .71 | 1.00 |
| 4 dyn unb | -.018 | +.009 | .69 | .75 | .74 | .75 | .74 | 1.00 |
Story: heteroskedasticity stresses the DYNAMIC cases (std SE undercovers,
SE/SD~0.71 -> robust/cluster help, BC reduces bias); static IFE is fairly robust
(honest). Factor selection unaffected (r.ok ~1). Output files tagged `_het`;
collate with `MC_ERR=hetero Rscript collate_mc.R` -> tab*_het.tex.

### How to run the complex MC on CERES
  rsync -av data-raw/hpc/{mc_xtife_hpc.R,ife.R,ife_unbalanced.R,collate_mc.R,submit_mc_xtife.sh} \
            bc25911@hpc.essex.ac.uk:/home/bc25911/xtife_mc/
  bash submit_mc_xtife.sh --hetero          # complex MC (writes *_het.rds)
  # download results, then:  MC_ERR=hetero Rscript collate_mc.R

---

## DGP4 (dynamic unbalanced) BC — KNOWN OPEN ISSUE: SWW-weak `b2` overcorrects under missingness

### Symptom (HPC full B=1000, hetero)
- d4 N=200 T=200 f80: IFE −0.0079 -> BC +0.0129 ; cov std .57 rob .74 cl .74 BC .49
- d4 N=200 T=100 f60: IFE −0.0151 -> BC +0.0150 ; cov std .53 rob .66 cl .67 BC .56
- d4 N=200 T=200 f60: IFE −0.0078 -> BC +0.0212 ; cov std .65 rob .79 cl .79 BC .31
BC flips the sign and worsens |bias|; coverage collapses. Note raw rob/cl cov ~.74-.79
is ALSO low — but that is the *uncorrected* Nickell bias (|bias|>SD at N=T=200), same
root cause: the dynamic bias is not being correctly removed.

### What was verified (so this is NOT any of these)
1. `b2` formula is FAITHFUL to SWW Thm 4.2(ii) (read clear PDF p.20): b2 corresponds to
   MW(2017) B1. b3+b5 = O_p(sqrt(T/N)) (Bai B / MW B2); b4+b6 = O_p(sqrt(N/T)) (Bai C / MW B3).
2. `b2` WORKS on BALANCED dynamic data: scales ~1/T, removes bias, matches the (fixed,
   working) MW B1 path. Decomp on balanced: b2 −1.04,−0.87,−0.72 for T=50,100,150 (decreasing ✓).
3. `b2` OVERCORRECTS only on UNBALANCED data: b2 flat/increasing (−1.42,−1.35,−1.82) ->
   correction decays slower than the 1/T raw bias -> overshoots, worst at T>=N + hetero.
4. b3+b5 cancel correctly even under hetero (sign fix holds); b4,b6 negligible. The whole
   overcorrection is `b2`.
5. `[L_lam_inv]_i` is correct: per-unit observed inverse == (1/Phi_i)(full inverse) in
   expectation; `T_i*|diag L_lam|`≈1 in both balanced and unbalanced. Replacing the per-unit
   inverse with the CLEAN full inverse × 1/Phi_i (paper-exact, = MW B1's shared Pf structure)
   changed b2 by <5% and did NOT fix the overcorrection -> bug is NOT the inverse/projection,
   and the MW-B1-structure rewrite is mathematically equivalent for MCAR (constant Phi_i).

### Conclusion
The overcorrection is a deeper FINITE-SAMPLE effect under missingness (most likely the
residual-bias feedback v_hat_it ⊃ −x_it(beta_hat−beta), amplified by the larger raw bias
when data are dropped, interacting with the kappa=sqrt(N/T) regime). It does NOT yield to
the inverse/projection rewrite. A definitive fix needs the SWW authors' reference code or a
dedicated finite-sample re-derivation — beyond inspection-based debugging.

### Status of the four BC paths
- DGP1 static balanced (Bai `.bias_correct`)            : FIXED, works.
- DGP2 static unbalanced (SWW strict `.ife_bias_unb`)   : works for N>=50 (near-no-op, small bias).
- DGP3 dynamic balanced (MW `.bias_correct_mw`)         : FIXED, works (cov ~.93).
- DGP4 dynamic unbalanced (SWW weak `b2`)               : OPEN — overcorrects under missingness.

### Recommendation for the JSS study
Present DGP1/2/3 as the BC results. For DGP4 report the RAW IFE (consistent; bias O(1/T)
shrinks) and flag BC as a documented limitation, OR restrict DGP4 to homoskedastic N>=T,
OR omit DGP4. Do NOT ship a guessed `b2` rescale (×kappa fixed one cell, broke others).

---

## VALIDATION against SWW Table 4 (replicate their exact DGP) — 2026-06-07

User's idea: replicate SWW's own DGP + Table 4 to prove the code. Implemented SWW DGP 2
(dynamic), Pattern 1 (MCAR p_it~U(0.5,0.9)): AR(1) factors rho_f=0.5 sig_f=0.5, lambda~N(1,1),
x1=y_{i,t-1}, x2~N(0,1), beta=(0.3,1), v_it~N(0,sig_v^2) sig_v~U(0.5,1.5), r=2, c=2 Bartlett.

Result (beta_1, the dynamic coef; B=300, se=hac exog=weak init=nnr):
              RAW bias   RAW sd   RAW CV |  BC bias   BC sd   BC CV
 OURS T=50:   -0.0088   0.0183   0.827  |  +0.0032  0.0219   0.803
 SWW  T=50:   -0.0084   0.0127   0.863  |  -0.0029  0.0125   0.924
 OURS T=100:  -0.0045   0.0116   0.883  |  +0.0037  0.0145   0.807
 SWW  T=100:  -0.0068   0.0085   0.857  |  -0.0021  0.0084   0.933

CONCLUSIONS:
1. RAW IFE ESTIMATOR IS CORRECT: raw bias matches SWW (-0.0088 vs -0.0084). VALIDATED.
2. BC b2 OVERCORRECTS ~2x: our correction +0.0120 vs SWW +0.0055 (T=50), ratio ~2.2x.
   Overshoots past 0 (BC bias +0.0032 vs SWW -0.0029) and inflates sd (x1.20 vs SWW x1.0).
3. Bandwidth (c=2), Bartlett kernel, b2 formula, [L_lam_inv] — ALL match the paper, yet 2x off.
   The 2x is NOT bandwidth/kernel/inverse. Likely a het/residual-feedback or pairing subtlety.
4. Larger RAW sd (0.0183 vs 0.0127) = NNR-vs-AM algorithm noise (Table 3 confirms NNR >> AM
   in RMSE); separate from the BC bug.
5. EVEN SWW's own BC is partial (never full) with CV 0.83-0.94 (<0.95) for dynamic unbalanced
   — so sub-nominal coverage is intrinsic to the method, not just our code.

Benchmark harness: /tmp/sww_repl.R pattern (SWW DGP2 P1). Use to hunt the 2x in b2.

---

## Appendix-G re-derivation OUTCOME (2026-06-07): b2 code is FAITHFUL to the paper

Read SWW p.20-21 as IMAGES (pdftotext was garbled) and verified b2 LINE-BY-LINE:
- b2k = (1/sqrt(NT)) sum_i sum_{t=1}^{T-1} sum_{s=t+1}^T Gamma((s-t)/L_T) d_is v_it d_it x_isk
        f_t'[L_lam_inv]_i f_s   (v at earlier t, x at later s, ONE-SIDED) -- code matches.
- [L_lam_inv]_i = (-sum_t d_it f_t f_t')^{-1} (minus inside) -- code matches.
- W_x = (1/NT) sum d xdot xdot'; beta_abc = beta - (1/sqrt(NT)) W_x^{-1} b_hat, b_hat=sum_{l=1}^6 b_l -- code matches.
- Kernel: theorem states TRUNCATION Gamma=1{|s-t|<=L_T}; sim text says BARTLETT. Code uses Bartlett
  = sim. (Their own paper formula vs sim text are inconsistent on the kernel.)
- L_T = floor(2 T^{1/5}) c=2 -- matches.

ORACLE TEST (decisive): computed b2 with TRUE factors + TRUE errors (projection f'[L]f = -[P_F] is
rotation/scale invariant, so true F usable). N=100 T=50 B=40:
  b2(code, estimated F,v) = -0.6458 ; b2(ORACLE true F,v) = -0.7123 ; ratio 1.10.
=> Fit quality is NOT the cause. The b2 FORMULA itself (perfect inputs) gives the ~2x correction.

CONCLUSION: our b2 is a CORRECT implementation of the published formula. The ~2x vs Table 4 is
NOT a bug in our code -- it reproduces the paper's formula exactly and over-corrects 2x even with
true inputs. The gap is therefore a SWW paper-vs-their-code discrepancy (kernel/damping/undocumented
factor) or an unobservable DGP detail. Not closable without SWW's reference code.

VALIDATED: estimator correct (raw bias matches Table 4); BC formula faithfully implemented.

---

## PUBLISHED version check (J.Econometrics 255 (2026) 106222) + full proofs — 2026-06-07

Read publishedunbalanced/ : main paper + supplement (mmc1, 34pp incl. full proofs).
NEW EVIDENCE (strengthens, does not change, the conclusion):
1. Published b2 formula == preprint == our code (verified main text p.11 images).
2. FULL PROOF in supplement Lemma E.3 derives the population bias b1+b2+b3+b4 (proof's b3=main b3+b5,
   proof's b4=main b4+b6). Proof's b2 = (1/sqrt NT) sum_i sum_t sum_{s=1}^T E(d_is v_is d_it x_it)
   f_t'[Lbar_lam']_i f_s -- sums ALL s, reduces to the one-sided main-text form (predetermined => nonzero
   only v-time<x-time; [Lbar] symmetric => weights match). So proof == main text == code. NO dropped factor.
3. Tested the dynamic missing ambiguity (Pattern 1): "pair" (d=dy_it*dy_{i,t-1}, Phi~0.49) vs "direct"
   (d~Bern(p_it), Phi~0.70, latent lag). BOTH overshoot vs Table 4 (BC +0.0010 / +0.0067 vs SWW -0.0029);
   direct is WORSE. So the missing mechanism does not explain the 2x.
4. No replication code published (only DOI 10.1016/j.jeconom.2026.106222 + supplement).

SUPERSEDED: the "2x over-correction" was caused by a bug found 2026-06-07 — see section below.

### A1b. The REAL root cause: `unique()` vs `sort(unique())` — BUG FOUND & FIXED (2026-06-07)

**Symptom.** The SWW b2 dynamic-bias term overcorrected by ~2× on unbalanced panels
(e.g., DGP4 at N=T=200: raw -0.008 → BC +0.013 instead of ~0). On balanced panels, b2
behaved correctly (matched Moon & Weidner B1).

**Root cause.** In `ife_unbalanced()` (line 1116), `time_vals` was computed as
`unique(data[[time_col]])`. In R, `unique()` preserves **first-appearance order**, not
sorted order. Since the data is sorted by (id, time), the first unique times come from
**unit 1's observed time periods** (which has gaps in unbalanced panels), then remaining
times from later units. For example, if unit 1 observes times {1,3,6,10,...}, then
`time_vals = {1,3,6,10,...,50, 2,4,5,...}` — NOT {1,2,3,...,50}.

This caused `time_idx = match(data$time, time_vals)` to produce **scrambled** indices.
When these indices were used in the b2 kernel computation (`gap = s_val - t_val`), the
"gap" was an index difference in the scrambled ordering, not the actual time gap. This
corrupted:
- **b2** (dynamic bias term, `exog = "weak"`) — ~15× inflation (not merely 2×; the
  previously-observed "2× overcorrection" was the net effect after partial cancellation)
- **HAC SE** (`se = "hac"`) — underestimation due to wrong kernel weights
- **HAC b4/b6** — similarly corrupted

Things NOT affected (because they don't use time gaps):
- The EM estimator and raw beta (order-invariant)
- i.i.d./robust/cluster SE
- b3, b4_iid, b5, b6_iid (use per-observation v²)
- delta/omega alternating LS (self-consistent under scrambled indexing)

**Fix.** One word: `sort()`.
```r
# Before (WRONG for unbalanced):
time_vals <- unique(data[[time_col]])
# After (CORRECT):
time_vals <- sort(unique(data[[time_col]]))
```
Applied in `ife_unbalanced()` and `ife_select_r_unb()` (both occurrences).

**Diagnosis.** Direct comparison of MW B1 (balanced, works) vs SWW b2 (unbalanced) on
the same balanced panel showed identical corrections (ratio 1.25, ≈1 within MCSE). On
70% unbalanced, SWW b2 exploded 15.6× — caused by 99/100 units having unsorted
`unit_obs_t`. After fix, all units have sorted `unit_obs_t` and b2 is well-behaved.

**Validation (40-rep local test, post-fix, DGP4 N=100 T=50 fill=80%).**
- Raw bias: -0.033, BC bias: -0.015 (55% reduction) ✅
- Raw coverage (std): 0.43 → BC+std: 0.71 → BC+HAC: 0.81
At T=100: Raw bias -0.018 → BC bias -0.007 (60% reduction), BC+HAC cov 0.84.

**Conclusion.** The "2× overcorrection" was NOT a discrepancy between SWW's formula and
their Table 4, nor a theoretical issue — it was a one-line ordering bug in our code.
The previous conclusion ("gap is internal to SWW") is RETRACTED.
`devtools::check()` passes: 0 errors, 0 warnings, 0 notes.
