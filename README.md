# Goal

Celina discovered an reproducibilty issue with clustering in her spatial data. In this repo, we track our investigations.

## Usage

Run this on terminal on the cluster:

```shell
bash sign.sh
```
## Backgound

Issue in `RunPCA` in Seurat when `approx = T`. Eventually calls `irlba` in R package `irlba`.
