#!/usr/bin/bash -l

module load rstudio

echo "Render script: " $1

# amd
module load amd
srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G \
  apptainer exec /apps/u24/opt/containers/rstudio/4.5.2.sif \
  quarto render $1 --output "amd.html"

# intel
module load intel
srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G \
  apptainer exec /apps/u24/opt/containers/rstudio/4.5.2.sif \
  quarto render $1 --output "intel.html"
