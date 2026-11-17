#!/bin/bash

# Run independent homologizer MCMC replicates in parallel on this machine
# (run id = seed), for convergence checks. Run from the project root:
#   ./scripts/homologizer/run_homologizer.sh [n_runs] [n_gen]    # default 2 runs, n_gen from homologizer.Rev
# Each replicate is single-threaded; screen output goes to log/homologizer_run<N>.out.

N_RUNS="${1:-2}"
N_GEN="${2:-}"    # optional: override n_gen (e.g. a short test run)

command -v rb >/dev/null || { echo "rb (RevBayes) not found on PATH"; exit 1; }
[ -f output/homologizer/rev/loci.Rev ] || { echo "run prep_homologizer.jl first"; exit 1; }

mkdir -p output/homologizer/mcmc log

pids=()
for i in $(seq 1 "${N_RUNS}"); do
  WRAPPER="output/homologizer/mcmc/run${i}_wrapper.Rev"
  { echo "run_id = ${i}"; [ -n "${N_GEN}" ] && echo "n_gen = ${N_GEN}"; echo 'source("scripts/homologizer/homologizer.Rev")'; } > "${WRAPPER}"
  echo "Run ${i} started $(date)"
  rb "${WRAPPER}" > "log/homologizer_run${i}.out" 2>&1 &
  pids+=($!)
done

status=0
for i in "${!pids[@]}"; do
  if wait "${pids[$i]}"; then
    echo "Run $((i + 1)) finished $(date)"
  else
    echo "Run $((i + 1)) FAILED, see log/homologizer_run$((i + 1)).out"; status=1
  fi
done
exit ${status}
