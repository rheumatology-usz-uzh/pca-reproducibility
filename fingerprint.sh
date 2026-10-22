#!/usr/bin/bash -l
# usage: bash fingerprint.sh amd      then: bash fingerprint.sh intel
# Runs blas_fingerprint.R natively and with forced kernels, at 1 and 2 OpenBLAS threads,
# inside one allocation on the requested node type; writes one log per condition.

NODE=${1:-amd}
SIF=/apps/u24/opt/containers/rstudio/4.5.2.sif

module load apptainer
module load rstudio
module load $NODE

srun -n 1 -c 2 --time=00:30:00 --mem=7G bash -c "
for CT in native Haswell Cooperlake; do
  for TH in 1 2; do
    ENV=OPENBLAS_VERBOSE=2,OPENBLAS_NUM_THREADS=\$TH
    [ \$CT != native ] && ENV=\$ENV,OPENBLAS_CORETYPE=\$CT
    apptainer exec --env \$ENV $SIF Rscript blas_fingerprint.R > fp_${NODE}_\${CT}_t\${TH}.txt 2>&1
  done
done"

# one line per run (all logs so far); identical fingerprints = bit-identical results
printf '%-26s %s\n' "run" "gemv         gemvT        gemm         syrk         svd          prcomp       irlba        iterations"
printf '%-26s %s\n' "sandbox_Cooperlake_t1" "051b6e92c983 e5fd2d526a13 1ecdbf089f5f c75fecd68edb 0bceb07aef9f 2dca42871ace 6c6519db7855 1d0b99180085"
printf '%-26s %s\n' "sandbox_Haswell_t1"    "7b75cba4c12b e5fd2d526a13 921f7d1ac5a8 ff7df8762430 cd2273e80c0e d4ac5af11100 656b49a67ef7 1d0b99180085"
for f in fp_*.txt; do
  printf '%-26s %s   %s\n' "${f%.txt}" \
    "$(grep -E '^(A |crossprod|svd|prcomp|irlba)' "$f" | awk '{print $NF}' | paste -sd' ')" \
    "$(grep -m1 '^Core:' "$f")"
done
