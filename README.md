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

## Next Steps

The principal components in PCA can have flipped signed because the sign is mathemtically not identifiable. The issue is that this happens even when fixing a seed. So far, we could reproduce by running it on AMD and Intel processors. Not reported on the issue page on `irlba` repo yet.

Within processor seems OK. Some todos:

- Check code in `irlba`
- Open issue on `irlba` repo
- Paper raising awarness of this kind of issues, use it as an example
