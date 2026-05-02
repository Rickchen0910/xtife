#!/bin/bash
# =============================================================================
# run_mc_job.sh  —  SGE job script for one (dgp_id, N, TT, fill) cell
#
# Arguments passed by submit_mc_xtife.sh:
#   $1  dgp_id  : 1/2/3/4
#   $2  N       : cross-sectional units
#   $3  TT      : time periods
#   $4  fill    : observed fraction ×100 (e.g. 80 for 80%)
#   $5  B       : replications (e.g. 1000)
# =============================================================================

#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -l mem_free=8G
#$ -m be
#$ -M bc25911@essex.ac.uk
#$ -pe smp 8

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
DGP=$1
NOBS=$2
TPER=$3
FILL=$4
BREPS=${5:-1000}

BASE=/home/bc25911/xtife_mc
RSCRIPT=$BASE/mc_xtife_hpc.R

mkdir -p "$BASE/results" "$BASE/logs/success" "$BASE/logs/failed"

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
echo "============================================================"
echo "  xtife Monte Carlo — DGP=$DGP  N=$NOBS  T=$TPER  fill=${FILL}%  B=$BREPS"
echo "  JOB_ID=$JOB_ID   JOB_NAME=$JOB_NAME"
echo "  Host: $(hostname)   Cores (NSLOTS): $NSLOTS"
echo "  Start: $(date)"
echo "============================================================"

# ---------------------------------------------------------------------------
# R setup
# ---------------------------------------------------------------------------
# xtife is loaded by mc_xtife_hpc.R by sourcing R/ife.R + R/ife_unbalanced.R
# from $BASE/xtife_src/ — no installation needed.
# (Uncomment the lines below only if you installed xtife to a user library
#  instead of using the source approach.)
# export R_LIBS_USER=/home/bc25911/R/library:$R_LIBS_USER

# Uncomment if CERES requires a module load for R:
# module load R/4.3.1

# ---------------------------------------------------------------------------
# Run R
# ---------------------------------------------------------------------------
/usr/bin/Rscript "$RSCRIPT" "$DGP" "$NOBS" "$TPER" "$FILL" "$NSLOTS" "$BREPS"
STATUS=$?

# ---------------------------------------------------------------------------
# Footer and log routing
# ---------------------------------------------------------------------------
echo "------------------------------------------------------------"
echo "  End: $(date)   Exit status: $STATUS"
echo "------------------------------------------------------------"

if [ $STATUS -eq 0 ]; then
    touch "$BASE/logs/success/${JOB_NAME}.${JOB_ID}"
else
    echo ">>> Job FAILED — moving log to logs/failed/"
    mv "$BASE/logs/out/${JOB_NAME}.out" "$BASE/logs/failed/" 2>/dev/null
fi

exit $STATUS
