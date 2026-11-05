#!/usr/bin/bash -l
# clustree_kernels.sh -- runs clustree_kernels.R on the ScienceCluster, from the project directory.
#
#   bash clustree_kernels.sh amd       # all kernel x thread configurations on one AMD node
#   bash clustree_kernels.sh intel     # the same on one Intel node
#   bash clustree_kernels.sh compare   # tables and clustering trees from all runs so far
#
# Override the allocation with CORES=128 MEM=150G bash clustree_kernels.sh amd
# One-time setup, to install clustree into the project's renv library:
#   apptainer exec /apps/u24/opt/containers/rstudio/4.5.2.sif Rscript -e 'renv::install("clustree"); renv::snapshot()'

NODE=${1:-amd}
CORES=${CORES:-256}
MEM=${MEM:-200G}
SIF=/apps/u24/opt/containers/rstudio/4.5.2.sif

module load apptainer
module load rstudio

# the tutorial data; download on the login node, since compute nodes may have no internet
if [ ! -d filtered_gene_bc_matrices/hg19 ]; then
  curl -L -o pbmc3k_filtered_gene_bc_matrices.tar.gz \
    https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
  tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
fi

if [ "$NODE" = compare ]; then
  srun -n 1 -c 4 --time=01:00:00 --mem=32G \
    apptainer exec $SIF Rscript clustree_kernels.R compare
else
  module load $NODE
  srun -n 1 -c $CORES --time=04:00:00 --mem=$MEM \
    apptainer exec $SIF Rscript clustree_kernels.R run --node=$NODE --cores=$CORES
fi
