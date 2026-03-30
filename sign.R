# testing sign flip issue in irlba

system("lscpu | grep 'Model name'")
system("ldd $(R RHOME)/lib/libR.so")

set.seed(1)

p <- 6988
n <- 477
x <- matrix(rnorm(p*n), nrow = n)

pca_prcomp <- prcomp(x = x, rank. = 30)
pca_irlba <- irlba::irlba(A = x, nv = 30)

head(pca_prcomp$rotation)
head(pca_irlba$v)

sessionInfo()
extSoftVersion()
