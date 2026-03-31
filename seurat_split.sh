#!/usr/bin/bash -l

module load rstudio

# amd
module load amd
srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G \
  apptainer exec /apps/u24/opt/containers/rstudio/4.5.2.sif \
  quarto render seurat_pca.qmd --output "amd_pca.html"

srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G \
  apptainer exec /apps/u24/opt/containers/rstudio/4.5.2.sif \
  quarto render seurat_cluster.qmd --output "amd_cluster.html"

# intel
module load intel
srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G \
  apptainer exec /apps/u24/opt/containers/rstudio/4.5.2.sif \
  quarto render seurat_cluster.qmd --output "intel_cluster.html"
