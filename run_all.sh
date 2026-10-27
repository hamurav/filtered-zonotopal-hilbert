#!/usr/bin/env bash
# Rerun the computations of the repository and record their output.
#
#   ./run_all.sh [--quick] [SECTION ...]
#
# SECTION is one of certificates, checks, tables, search, strata, cycles,
# generic_exp (default: all, in this order). Every command writes its output
# to the log file named in README.md, next to the script, and one summary
# line (status, seconds, command) to $RUN_LOG (default run_all.log).
# --quick skips the runs marked "long" below (about 3 hours in total on
# four cores); the remaining runs take about 25 minutes.
# JOBS (default 4) is the number of worker processes for the scripts that
# accept --jobs.
#
# Requirements: Python 3.9 or later with NumPy, SymPy and networkx
# (requirements.txt); a C++17 compiler and GMP for certificates/.

set -u
cd "$(dirname "$0")"
ROOT=$(pwd)
QUICK=0
SECTIONS=()
for a in "$@"; do
  case "$a" in
    --quick) QUICK=1 ;;
    certificates|checks|tables|search|strata|cycles|generic_exp) SECTIONS+=("$a") ;;
    *) echo "unknown argument: $a" >&2; exit 2 ;;
  esac
done
[ ${#SECTIONS[@]} -eq 0 ] && SECTIONS=(certificates checks tables search strata cycles generic_exp)
JOBS=${JOBS:-4}
RUN_LOG=${RUN_LOG:-$ROOT/run_all.log}
FAILED=0

# run [long] DIR LOG COMMAND...  (LOG is relative to DIR; '-' if the script
# writes its own output and only its stdout is discarded)
run() {
  local long=0
  if [ "$1" = long ]; then long=1; shift; fi
  local dir=$1 log=$2; shift 2
  if [ $QUICK = 1 ] && [ $long = 1 ]; then
    printf 'skipped        %s: %s\n' "$dir" "$*" | tee -a "$RUN_LOG"
    return
  fi
  local out=/dev/null
  [ "$log" != - ] && out="$ROOT/$dir/$log"
  local t0 t1 st
  t0=$(date +%s)
  if (cd "$ROOT/$dir" && "$@") > "$out" 2>&1; then st=ok; else st=FAILED; FAILED=1; fi
  t1=$(date +%s)
  printf '%-6s %7ss  %s: %s%s\n' "$st" $((t1 - t0)) "$dir" "$*" \
    "$([ "$log" != - ] && echo "  > $log")" | tee -a "$RUN_LOG"
}

py() { OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 "$@"; }

echo "# $(date -u '+%Y-%m-%d %H:%M UTC')  $(uname -srm)  $(python3 --version 2>&1)  sections: ${SECTIONS[*]}$([ $QUICK = 1 ] && echo ' (quick)')" | tee -a "$RUN_LOG"

for s in "${SECTIONS[@]}"; do
  case $s in
  certificates)
    run certificates verified_results.txt py verify_parameter_results.py
    run certificates tree_messages.log py tree_messages.py
    ;;
  checks)
    run checks t3_profile.log py t3_profile.py ../certificates ../strata/d4_witness_lib.py
    run checks exp_check.log py exp_check.py ../certificates
    run checks d6e6_check.log py d6e6_check.py ../certificates
    run checks c6_annihilator.log py c6_annihilator.py
    run checks neck_states.log py neck_states.py ../certificates
    run checks neck_theta8.log py neck_direct.py theta 8
    run checks neck_K13.log py neck_direct.py K13 6
    run checks neck_K3.log py neck_direct.py K3 6
    run checks paths_tails.log py paths_tails.py
    run checks loci_algebra.log py loci_algebra.py
    run long checks subdivision.log py subdivision.py
    ;;
  tables)
    if [ $QUICK = 1 ]; then
      run tables check_tables_mod.txt py check_tables.py tables345_clean.tex --certificates ../certificates --mod 1048573
    else
      run tables check_tables_output.txt py check_tables.py tables345_clean.tex --certificates ../certificates --json tables.json
    fi
    run long tables complete_profiles.log py complete_profiles.py tables345_clean.tex
    ;;
  search)
    run search results/small_graphs.log py small_graphs.py
    run search results/pairs_run.log py verify_pairs.py --jobs "$JOBS"
    if [ $QUICK = 1 ]; then
      run search results/trees_quick_run.log py search.py trees --nmax 10 --jobs "$JOBS" --tag _quick
    else
      run search results/trees_run.log py search.py trees --nmax 10 --direct --jobs "$JOBS"
    fi
    run search results/trees_without_u4_run.log py search.py trees --nmax 10 --jobs "$JOBS" --tag _without_u4 --functions u u+u^2 u+u^2+u^3 u+u^2+2u^3 u+u^3
    run search results/unicycles_run.log py search.py unicycles --nmax 8 --jobs "$JOBS"
    run search results/nsh_example2.txt py nsh_example2.py
    run search results/exp_small.log py exp_small.py --jobs "$JOBS"
    run search - py degenerate_check.py trees --nmax 12
    run search - py degenerate_check.py lemma
    run search - py degenerate_check.py cycles
    ;;
  strata)
    run strata d3_closure.log py d4_closure.py ../certificates 3
    run strata d4_closure.log py d4_closure.py ../certificates 4
    run strata d4_cells.log py d4_cells.py ../certificates
    run strata d4_witness.log py d4_witness.py ../certificates 10
    run long strata d5_closure.log py d5_closure.py ../certificates
    run long strata surf_generic.log py d5_surfaces.py ../certificates generic
    for i in 0 1 2 3 4 5 6 7 8; do
      run long strata surf_$i.log py d5_surfaces.py ../certificates $i
    done
    run long strata d5_witness.log py d5_witness.py ../certificates 14
    run long strata collisions_D3.log py collisions.py ../certificates 3 18
    run long strata collisions_D4.log py collisions.py ../certificates 4 14
    run long strata info_D3.log py info_per_point.py ../certificates 3 18
    run long strata info_D4.log py info_per_point.py ../certificates 4 14
    run long strata moves_D3.log py collision_moves.py ../certificates 3 16
    run long strata swaps_D3.log py swap_test.py ../certificates 3 14
    ;;
  cycles)
    if [ $QUICK = 1 ]; then
      run cycles results/cycle_hilbert_quick.log python3 cycle_hilbert.py 3 13 --primes 2 --out results/cycle_series_quick.json
      run cycles results/cycle_formulas_quick.log py cycle_formulas.py results/cycle_series_quick.json
    else
      run cycles results/cycle_hilbert.log python3 cycle_hilbert.py 3 17 --primes 2 --out results/cycle_series.json
      run cycles results/cycle_formulas.log py cycle_formulas.py results/cycle_series.json
    fi
    run cycles results/annihilators.log py annihilators.py
    ;;
  generic_exp)
    run checks/generic_exp results/phase1.log py phase1.py
    run long checks/generic_exp results/scan3_p101.log py scan.py 101 ALL3
    run long checks/generic_exp results/scan4_p101.log py scan.py 101 ALL4
    run long checks/generic_exp results/scan_p61.log py scan.py 61 ALL3 ALL4
    run long checks/generic_exp results/classify_p101.txt py classify.py 101
    run long checks/generic_exp results/classify_p61.txt py classify.py 61
    run long checks/generic_exp results/six_6.log py six.py 6 15
    run checks/generic_exp results/bonus_low.log py bonus_low.py
    run checks/generic_exp - py low_certificate.py --jobs "$JOBS"
    run long checks/generic_exp results/summary.txt py summary.py
    ;;
  esac
done

if [ $FAILED = 0 ]; then
  echo "# all runs finished without error" | tee -a "$RUN_LOG"
else
  echo "# SOME RUNS FAILED (see the lines marked FAILED)" | tee -a "$RUN_LOG"
fi
exit $FAILED
