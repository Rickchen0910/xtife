#!/bin/bash
# =============================================================================
# submit_mc_xtife.sh
# Submits one SGE job per (dgp_id, N, TT, fill) parameter cell.
#
# Usage:
#   bash submit_mc_xtife.sh          # full run, B = 1000 reps
#   bash submit_mc_xtife.sh --test   # smoke test, B = 10 reps
#
# Parameter grid
# -----------------------------------------------------------------------
#  dgp_id  : 1 = static balanced   (ife, Bai 2009 BC)
#             2 = static unbalanced (ife_unbalanced, SWW2025 strict BC)
#             3 = dynamic balanced  (ife, Moon & Weidner 2017 BC)
#             4 = dynamic unbalanced(ife_unbalanced, SWW2025 weak BC)
#  N       : 30 50 100 200
#  TT      : 10 50 100 200  (T=10 skipped for dynamic DGPs 3 & 4)
#  fill    : 80 60           (only varies for unbalanced DGPs 2 & 4;
#                             balanced DGPs submit once with fill=80)
# =============================================================================

BASE=/home/bc25911/xtife_mc
JOB_SCRIPT=$BASE/run_mc_job.sh

mkdir -p "$BASE/logs/out" "$BASE/logs/err" \
         "$BASE/logs/success" "$BASE/logs/failed" \
         "$BASE/results" "$BASE/xtife_src"

# xtife is loaded by sourcing R files directly — no package installation needed.
# Before running this script, upload the two source files once:
#   rsync -av R/ife.R R/ife_unbalanced.R \
#             bc25911@hpc.essex.ac.uk:$BASE/xtife_src/
if [ ! -f "$BASE/xtife_src/ife.R" ]; then
  echo "WARNING: $BASE/xtife_src/ife.R not found."
  echo "  Upload with: rsync -av R/ife.R R/ife_unbalanced.R \\"
  echo "    bc25911@hpc.essex.ac.uk:$BASE/xtife_src/"
  echo "  Continuing submission — jobs will fail if files are missing."
fi

# ---------------------------------------------------------------------------
# Flags (order-independent):
#   --test       -> B=10 smoke test
#   --hetero     -> heteroskedastic errors          (MC_ERR=hetero)
#   --serial     -> MA(1) serial correlation        (MC_ERR=serial)
#   --serialhet  -> serial correlation + heterosk.  (MC_ERR=serialhet)
# (serial correlation applies to the static DGPs 1 & 2 only; the dynamic DGPs
#  stay serially uncorrelated — see mc_xtife_hpc.R.)
# ---------------------------------------------------------------------------
B=1000
MCERR=homo
for arg in "$@"; do
  case "$arg" in
    --test)      B=10 ;;
    --hetero)    MCERR=hetero ;;
    --serial)    MCERR=serial ;;
    --serialhet) MCERR=serialhet ;;
  esac
done
if [[ $B -eq 10 ]]; then echo ">>> TEST MODE (B=10)"; else echo ">>> FULL RUN (B=1000)"; fi
echo ">>> Error structure: MC_ERR=$MCERR"
# job-name / log tag so different error structures do not collide
case "$MCERR" in
  hetero)    JT="_het" ;;
  serial)    JT="_ser" ;;
  serialhet) JT="_serhet" ;;
  *)         JT="" ;;
esac

# ---------------------------------------------------------------------------
# Parameter arrays
# ---------------------------------------------------------------------------
# DGPs to run: override with the DGPS env var, e.g.
#   DGPS="2 4" bash submit_mc_xtife.sh        # only the unbalanced DGPs
dgps=(${DGPS:-1 2 3 4})
obs=(30 50 100 200)
tpers=(10 50 100 200)
fills=(80 60)          # ×100; only relevant for DGPs 2 & 4

# ---------------------------------------------------------------------------
# OLS vs NNR comparison experiment (Section 6)
# Fixed: N=50 T=50 fill=50 r_true=3 — submit as special dgp=5
# ---------------------------------------------------------------------------
# Uncomment the block below once mc_xtife_hpc.R supports dgp_id=5
# qsub -N mc_xtife_d5_N50_T50_f50 \
#      -o "$BASE/logs/out/mc_xtife_d5_N50_T50_f50.out" \
#      -e "$BASE/logs/err/mc_xtife_d5_N50_T50_f50.err" \
#      "$JOB_SCRIPT" 5 50 50 50 "$B"

# ---------------------------------------------------------------------------
# Submission loop
# ---------------------------------------------------------------------------
n_submitted=0

for dgp in "${dgps[@]}"; do
  for nobs in "${obs[@]}"; do
    for tper in "${tpers[@]}"; do

      # Dynamic DGPs (3 & 4): burn-in requires T >= 50
      if [[ ( $dgp -eq 3 || $dgp -eq 4 ) && $tper -eq 10 ]]; then
        continue
      fi

      for fill in "${fills[@]}"; do

        # Balanced DGPs (1 & 3): fill rate is irrelevant — submit once only
        if [[ ( $dgp -eq 1 || $dgp -eq 3 ) && $fill -ne 80 ]]; then
          continue
        fi

        jobname="mc_xtife_d${dgp}_N${nobs}_T${tper}_f${fill}${JT}"

        qsub \
          -N  "$jobname" \
          -o  "$BASE/logs/out/${jobname}.out" \
          -e  "$BASE/logs/err/${jobname}.err" \
          -v  MC_ERR="$MCERR" \
          "$JOB_SCRIPT" \
          "$dgp" "$nobs" "$tper" "$fill" "$B"

        echo "  Submitted: $jobname"
        n_submitted=$(( n_submitted + 1 ))

      done   # fill
    done     # tper
  done       # nobs
done         # dgp

echo "---"
echo "Total jobs submitted: $n_submitted"
echo "Monitor with:  qstat -u bc25911"
echo "Results will appear in: $BASE/results/"
