#!/usr/bin/bash -l

module load rstudio

# amd
module load amd
srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G Rscript sign.R

# intel
module load intel
srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G Rscript sign.R
