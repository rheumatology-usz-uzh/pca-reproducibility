# make_figures.R -- reproduces Figure 1 (figs/agreement.pdf), the numbers in Table 2 and the
# summary statistics quoted in the "Numerical illustration" section of paper.tex.
# Run from the directory containing paper.tex:  Rscript make_figures.R
library(ggplot2)

set.seed(2026)
reps <- 40000
a0 <- 1.5 * 2^20                   # running total in the binade [2^20, 2^21)
h <- 2^(20 - 52)                   # its grid spacing: one ulp
ns <- c(1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000)
rterm <- function(k) (runif(k) + runif(k) * 2^-32) * 2^-10   # full-precision terms

# Theorem (two kernels): A adds each term to the large total, B sums at small magnitude first
sim <- do.call(rbind, lapply(ns, function(n) {
  sA <- rep(a0, reps)
  tB <- rep(0, reps)
  for (j in seq_len(n)) {
    x <- rterm(reps)
    sA <- sA + x
    tB <- tB + x
  }
  D <- (sA - (a0 + tB)) / h        # difference in ulps, always an integer
  stopifnot(all(D == round(D)))
  data.frame(n = n, p = mean(D == 0), se = sd(D == 0) / sqrt(reps), ED2 = mean(D^2))
}))

# Proposition (same grid): forward vs reverse summation onto the large total
same_grid <- do.call(rbind, lapply(c(10, 100, 1000), function(n) {
  X <- matrix(rterm(5000 * n), 5000, n)
  fwd <- rep(a0, 5000); for (j in seq_len(n)) fwd <- fwd + X[, j]
  bwd <- rep(a0, 5000); for (j in rev(seq_len(n))) bwd <- bwd + X[, j]
  data.frame(n = n, agree = mean(fwd == bwd))
}))

# f_k(0): density at 0 of the sum of k uniforms on (-1/2, 1/2), i.e. the cardinal B-spline of
# order k at its centre, computed with the Cox-de Boor recurrence (all terms non-negative)
f0 <- function(k) {
  x <- k / 2 - (0:(k - 1))
  v <- as.numeric(x >= 0 & x < 1)
  if (k == 1) return(v[1])
  for (j in 2:k) {
    i <- seq_len(k - j + 1)
    v <- (x[i] * v[i] + (j - x[i]) * v[i + 1]) / (j - 1)
  }
  v[1]
}

# Figure 1
kk <- unique(round(10^seq(log10(2), log10(2001), length.out = 60)))
theory <- rbind(
  data.frame(k = kk, p = sapply(kk, f0), curve = "Exact"),
  data.frame(k = kk, p = sqrt(6 / (pi * kk)), curve = "Diaconis-Freedman approximation")
)
fig <- ggplot() +
  geom_line(data = theory, aes(k, p, linetype = curve)) +
  geom_pointrange(data = sim, aes(n + 1, p, ymin = p - 2 * se, ymax = p + 2 * se,
                                  shape = "Simulation"), size = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  scale_linetype_manual(values = c("Exact" = "solid", "Diaconis-Freedman approximation" = "dashed")) +
  scale_shape_manual(values = c("Simulation" = 16)) +
  labs(x = "n + 1", y = "P(kernels agree bit for bit)", linetype = NULL, shape = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
dir.create("figs", showWarnings = FALSE)
ggsave("figs/agreement.pdf", fig, width = 6.5, height = 4)

# Numbers quoted in the text and in Table 2
cat(sprintf("max deviation from f_{n+1}(0): %.1f standard errors\n",
            max(abs(sim$p - sapply(sim$n + 1, f0)) / pmax(sim$se, 1e-12))))
cat("same-grid agreement (n = 10, 100, 1000):", paste0(round(100 * same_grid$agree, 2), "%"), "\n")
m <- subset(sim, n >= 2)
print(data.frame(n = m$n, simulated = round(m$ED2, 2), theory = round((m$n + 1) / 12, 2)),
      row.names = FALSE)
