#!/usr/bin/bash -l

module load apptainer
module load rstudio
module load amd

srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G \
	apptainer exec --env OPENBLAS_CORETYPE=Haswell,OPENBLAS_VERBOSE=2 \
        /apps/u24/opt/containers/rstudio/4.5.2.sif \
	Rscript demo.R

srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G \
        apptainer exec --env OPENBLAS_CORETYPE=Cooperlake,OPENBLAS_VERBOSE=2 \
        /apps/u24/opt/containers/rstudio/4.5.2.sif \
        Rscript demo.R

srun --pty -n 1 -c 2 --time=01:00:00 --mem=7G \
        apptainer exec --env OPENBLAS_VERBOSE=2 \
        /apps/u24/opt/containers/rstudio/4.5.2.sif \
        Rscript demo.R
