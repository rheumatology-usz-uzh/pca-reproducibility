# illustration.R -- reproduces the agreement figure (figs/agreement.pdf) and the numbers quoted
# in the "Numerical illustration" section of blas_kernels.tex.
# Run from the directory containing blas_kernels.tex:  Rscript illustration.R
library(ggplot2)

set.seed(2026)
reps <- 40000
a0 <- 1.5 * 2^20                   # running total in the binade [2^20, 2^21)
h <- 2^(20 - 52)                   # its grid spacing: one ulp
ns <- c(1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000)
rterm <- function(k) (runif(k) + runif(k) * 2^-32) * 2^-10   # full-precision terms

# Theorem 1: kernel A adds each term to the large running total,
# kernel B sums the terms at small magnitude first and adds the result once
sim <- do.call(rbind, lapply(ns, function(n) {
  sA <- rep(a0, reps)
  tB <- rep(0, reps)
  for (j in seq_len(n)) {
    x <- rterm(reps)
    sA <- sA + x
    tB <- tB + x
  }
  D <- (sA - (a0 + tB)) / h        # difference in ulps
  stopifnot(all(D == round(D)))    # Theorem 1 (i): always a whole number of ulps
  data.frame(n = n, p = mean(D == 0), se = sd(D == 0) / sqrt(reps))
}))
sim$approx <- sqrt(6 / (pi * sim$n))           # normal approximation, Theorem 1 (ii)

# Within one binade the order does not matter: forward vs reverse summation onto the large total
same_grid <- do.call(rbind, lapply(c(10, 100, 1000), function(n) {
  X <- matrix(rterm(5000 * n), 5000, n)
  fwd <- rep(a0, 5000); for (j in seq_len(n)) fwd <- fwd + X[, j]
  bwd <- rep(a0, 5000); for (j in rev(seq_len(n))) bwd <- bwd + X[, j]
  data.frame(n = n, agree = mean(fwd == bwd))
}))

# Figure
nn <- exp(seq(log(2), log(2000), length.out = 100))   # the approximation exceeds 1 at n = 1
fig <- ggplot() +
  geom_line(data = data.frame(n = nn, p = sqrt(6 / (pi * nn))),
            aes(n, p, linetype = "Normal approximation")) +
  geom_pointrange(data = sim, aes(n, p, ymin = p - 2 * se, ymax = p + 2 * se,
                                  shape = "Simulation"), size = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  scale_linetype_manual(values = c("Normal approximation" = "dashed")) +
  scale_shape_manual(values = c("Simulation" = 16)) +
  labs(x = "n (number of terms)", y = "P(kernels agree bit for bit)", linetype = NULL, shape = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
dir.create("figs", showWarnings = FALSE)
ggsave("figs/agreement.pdf", fig, width = 6.5, height = 4)

# Numbers quoted in the text
two <- sim$n == 2; ten <- sim$n == 10
cat(sprintf("n = 2: agreement rate %.3f (standard error %.3f); exact value 3/4; approximation %.3f\n",
            sim$p[two], sim$se[two], sim$approx[two]))
cat(sprintf("n = 10: the approximation exceeds the simulated rate by %.0f%%\n",
            100 * (sim$approx[ten] / sim$p[ten] - 1)))
z <- (sim$approx - sim$p) / sim$se
cat(sprintf("n >= 50: at most %.1f standard errors between approximation and simulation\n",
            max(abs(z[sim$n >= 50]))))
cat(sprintf("Discussion, n = 10: Diaconis-Freedman sqrt(6/(pi (n+1))) = %.3f, simulated %.3f, sqrt(6/(pi n)) = %.3f\n",
            sqrt(6 / (pi * 11)), sim$p[ten], sim$approx[ten]))
cat("same-grid agreement (n = 10, 100, 1000):", paste0(round(100 * same_grid$agree, 2), "%"), "\n")
