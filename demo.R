#!/usr/bin/env Rscript

set.seed(1)

A <- matrix(rnorm(2000 * 2000), 2000)  # rnorm() doesn't use BLAS, so A and v
v <- rnorm(2000)                        # are bit-identical in every run
x <- A %*% v                            # this product is computed by OpenBLA

print(x[1], digits = 22)

cpu <- try(system("lscpu | grep 'Model name'", intern = TRUE))
cpu
