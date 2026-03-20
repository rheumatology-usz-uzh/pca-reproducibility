# Goal

Celina discovered an reproducibilty issue with clustering in her spatial data. In this repo, we track our investigations.

## Usage

Run this on terminal on the cluster:

```shell
bash sign.sh
```
## Backgound

Issue in `RunPCA` in Seurat when `approx = T`. Eventually calls `irlba` in R package `irlba`.

This happens in `R/dimensional_reduction.R` [here](https://github.com/satijalab/seurat/blob/3f6462bbdc4e3d78e5c070ed0fc512b5a5bc3351/R/dimensional_reduction.R#L944-L945).
