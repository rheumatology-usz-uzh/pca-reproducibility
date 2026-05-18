#!/usr/bin/bash -l
#SBATCH --job-name=pca
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=7G

# input parameter: 'amd' or 'intel'
CPU=$1
echo "Load processor:" $CPU
module load $CPU
module load rstudio

Rscript -e "quarto::quarto_render( \
  input = 'seurat.qmd', \
  output_file = '${CPU}_exact.html', \
  execute_params = list(pca = 'exact') \
)"

Rscript -e "quarto::quarto_render( \
  input = 'seurat.qmd', \
  output_file = '${CPU}_irlba.html', \
  execute_params = list(pca = 'irlba') \
)"

Rscript -e "quarto::quarto_render( \
  input = 'seurat.qmd', \
  output_file = '${CPU}_rsvd.html', \
  execute_params = list(pca = 'rsvd') \
)"
