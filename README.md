# Reproducibility Issue in Seurat

Celina discovered an reproducibilty issue with clustering in her spatial data. In this repo, we track our investigations.

The principal components in PCA can have flipped signed because the sign is mathematically not identifiable. The issue is that this happens even when fixing a seed. So far, we could reproduce by running it on AMD and Intel processors. Not reported on the issue page on `irlba` repo yet. Within processor seems OK. 

## Usage

Run this on terminal on the cluster:

```shell
bash seurat_comparison.sh
```

## Backgound

Issue in `RunPCA` in Seurat when `approx = T`. Eventually calls `irlba` in R package `irlba`.

This happens in `R/dimensional_reduction.R` [here](https://github.com/satijalab/seurat/blob/3f6462bbdc4e3d78e5c070ed0fc512b5a5bc3351/R/dimensional_reduction.R#L944-L945).

## Related Work

- Found reproducibility issue in scanpy for [UMAP and Leiden analysis](https://github.com/grst/scanpy_reproducibility)
- Also mentioned in PCA [tutorial](https://ehrapy.readthedocs.io/en/latest/tutorials/notebooks/mimic_2_introduction.html#principle-component-analysis-pca) in ehrapy using scanpy

## Workaround

Current solution is to set `approx = F`.

## Next Steps

- [x] Check code in `irlba` -> not a bug, looks like it's a low-level difference in the linear algebra libraries BLAS and LAPACK
- [x] Reproduce issue in Seurat clustering [tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) by flipping sign after PCA, then running `FindNeighbors` and `FindClusters`. The clusters should be invariant to sign flips in PCs. -> Flipping signs doesn't affect UMAP or clustering. The issue must be in the numerical errors in `irlba`, not the signs.
- [x] Find paper on numerical error propagation in the iterative algorithm used in `irlba`. -> Didn't find anything.
- [x] Reproduce issue in data integration from OSCA [book](https://bioconductor.org/books/release/OSCA.multisample/integrating-datasets.html#mnn-correction)
- [x] Send abstract to conference on [Symposium on Meta Science for Methods Research](https://crsuzh.pages.uzh.ch/msmr/)
- [x] Replace openBLAS and openLAPACK libraries with Netlib solves the issue, but will be much solver
- [ ] Write paper to raise awareness of this kind of issues:
    - Motivate using Seurat and `approx = T` as default setting.
    - Give numerical illustration of the propagation of error.
    - Use ideas from statistical literature on round-off errors, e.g., [Diaconis and Freedman](https://drive.google.com/file/d/1GjD5WARrmAuF4wN0IJGEuCBK9judpcD4/view).
