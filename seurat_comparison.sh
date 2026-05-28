#!/usr/bin/bash -l

module load apptainer
module load rstudio

BLAS="/usr/lib/x86_64-linux-gnu/blas/libblas.so"
LAPACK="/usr/lib/x86_64-linux-gnu/lapack/liblapack.so"

CPU_TYPES=("amd" "intel")
PCA_ALGOS=("exact" "irlba" "rsvd")

for CPU in ${CPU_TYPES[@]}; do

    echo "load processor:" $CPU
    module load $CPU

    for PCA in ${PCA_ALGOS[@]}; do

      echo "pca algorithm:" $PCA
      srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G \
        apptainer exec --env LD_PRELOAD=$BLAS:$LAPACK \
        /apps/u24/opt/containers/rstudio/4.5.2.sif \
        quarto render seurat.qmd --output ${CPU}_${PCA}.html -P pca:$PCA

    done

done

echo "comparison"
srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G \
  apptainer exec /apps/u24/opt/containers/rstudio/4.5.2.sif \
  quarto render seurat_comparison.qmd

