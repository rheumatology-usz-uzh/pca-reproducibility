#!/usr/bin/bash -l

module load miniforge3
source activate general

# amd
module load amd
srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G \
  papermill pca_scikit.ipynb pca_scikit_amd.ipynb

# intel
module load intel
srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G \
  papermill pca_scikit.ipynb pca_scikit_intel.ipynb
