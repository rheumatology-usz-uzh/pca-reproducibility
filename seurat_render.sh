#!/usr/bin/bash -l

module load rstudio

CPU_TYPES=("amd" "intel")
PCA_ALGOS=("exact" "irlba" "rsvd")

for CPU in ${CPU_TYPES[@]}; do

    echo "load processor:" $CPU
    module load $CPU

    for PCA in ${PCA_ALGOS[@]}; do

        echo "pca algorithm:" $PCA

        srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G \
          Rscript -e "quarto::quarto_render( \
            input = 'seurat.qmd', \
            output_file = '${CPU}_${PCA}.html', \
            execute_params = list(pca = '${PCA}') \
          )"

    done

done

echo "comparison"
srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G \
  Rscript -e "quarto::quarto_render( \
    input = 'seurat_comparison.qmd') \
  )"
