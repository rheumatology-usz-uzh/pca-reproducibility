# blas_fingerprint.R -- which OpenBLAS kernel, CPU and thread setting was used, plus an
# md5 fingerprint of several BLAS/LAPACK results on fixed synthetic data.
# Same fingerprint on two runs = bit-identical results for that operation.
#
#   OPENBLAS_VERBOSE=2 OPENBLAS_NUM_THREADS=2 Rscript blas_fingerprint.R
#   OPENBLAS_VERBOSE=2 OPENBLAS_NUM_THREADS=2 OPENBLAS_CORETYPE=Haswell Rscript blas_fingerprint.R
#   inside apptainer: apptainer exec --env OPENBLAS_CORETYPE=Haswell,OPENBLAS_NUM_THREADS=2 ...

fp <- function(x) {                                   # md5 of the raw bytes, first 12 hex digits
  f <- tempfile(); on.exit(unlink(f)); writeBin(as.vector(unlist(x)), f)
  substr(unname(tools::md5sum(f)), 1, 12)
}
cpu <- sub(".*:\\s*", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1])
cat("cpu:       ", cpu, "\n")
ct <- Sys.getenv("OPENBLAS_CORETYPE"); cat("coretype:  ", if (ct == "") "(auto)" else ct, "\n")
cat("threads:   ", Sys.getenv("OPENBLAS_NUM_THREADS", "(unset)"),
    "| cpus in affinity mask:", length(parallel::mcaffinity()), "\n")
cat("blas:      ", normalizePath(sessionInfo()$BLAS), "\n")
cat("lapack:    ", normalizePath(sessionInfo()$LAPACK), "\n\n")

set.seed(1)                                           # rnorm() does not use BLAS: same inputs everywhere
A <- matrix(rnorm(2000 * 2000), 2000); v <- rnorm(2000)
B <- matrix(rnorm(2000 * 50), 2000)
S <- matrix(rnorm(600 * 400), 600)

res <- list(
  "A %*% v          (gemv)"   = A %*% v,
  "crossprod(A, v)  (gemv T)" = crossprod(A, v),
  "A %*% B          (gemm)"   = A %*% B,
  "crossprod(S)     (syrk)"   = crossprod(S),
  "svd(S)           (LAPACK)" = svd(S),
  "prcomp(S)"                 = prcomp(S)$x
)
if (requireNamespace("irlba", quietly = TRUE)) {
  set.seed(23); r <- irlba::irlba(S, nv = 20)
  res[["irlba(S, nv = 20)"]] <- list(r$d, r$u, r$v)
  res[["irlba iterations"]]  <- c(r$iter, r$mprod)
}
for (nm in names(res)) cat(sprintf("%-27s %s\n", nm, fp(res[[nm]])))
cat(sprintf("\nx[1] = %.17g\n", res[[1]][1]))
