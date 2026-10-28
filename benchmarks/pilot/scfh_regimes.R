# Pilot of scfh_sae() on the design of Maranzano, Mattera and Sugasawa (2026),
# sec. 3.1, Algorithm 1, Scenario 1 (bands along one coordinate, clear
# separation), on a 16 x 16 lattice with rook neighbours instead of the 256
# Po Valley areas (technical-report/DISCREPANCIES.md, SCFH-3).
#
# K = 3, beta0 = (50, 75, 100), beta1 = (-5, 2, 10), x ~ N(0, 4); R2 shares
# 0.90 / 0.07 / 0.03 (fixed / random / sampling). The paper does not give the
# formula for the variances from the R2 shares; here, per cluster,
# sigma2_uk = beta1k^2 * 4 * 0.07 / 0.90 and psi_d = beta1k^2 * 4 * 0.03 / 0.90.
#
# Part 1: ARI, share of correct labels and RMSE relative to the pooled FH
#   (K = 1, adjusted REML) at the true K = 3, by phi (paper: Table S2).
# Part 2: bootstrap MSE against the Monte Carlo MSE at K = 3, phi = 0.5
#   (paper: Table S16 for Scenarios 3-4).
#
# Uses the installed package. Rscript benchmarks/pilot/scfh_regimes.R [R1] [R2] [B] [cores]

suppressPackageStartupMessages(library(flexsae))
args <- commandArgs(trailingOnly = TRUE)
R1 <- as.integer(args[1] %||% 200)
R2 <- as.integer(args[2] %||% 100)
B <- as.integer(args[3] %||% 100)
cores <- as.integer(args[4] %||% 2)
RNGkind("L'Ecuyer-CMRG")
set.seed(20261028)

side <- 16
g <- expand.grid(col = seq_len(side), row = seq_len(side))
W <- 1 * (as.matrix(dist(g)) == 1)
regime <- cut(g$col, c(0, 5, 11, 16), labels = FALSE)
b0 <- c(50, 75, 100)
b1 <- c(-5, 2, 10)
v_fe <- b1^2 * 4
s2u <- v_fe * 0.07 / 0.90
psi <- v_fe[regime] * 0.03 / 0.90

draw <- function() {
  x <- rnorm(nrow(g), sd = 2)
  theta <- b0[regime] + b1[regime] * x + rnorm(nrow(g), sd = sqrt(s2u[regime]))
  data.frame(x = x, vardir = psi, theta = theta, y = theta + rnorm(nrow(g), sd = sqrt(psi)))
}
share_correct <- function(cl) {
  map <- vapply(1:3, function(j) {
    i <- which(cl == j)
    if (length(i) == 0) NA_integer_ else which.max(tabulate(regime[i], 3))
  }, integer(1))
  mean(map[cl] == regime)
}
ari <- flexsae:::.adjusted_rand
rmse <- function(fit, d) sqrt(mean((fit$df_est$est - d$theta)^2))

# Part 1
phis <- c(0, 0.25, 0.5, 1, 2)
part1 <- parallel::mclapply(seq_len(R1), function(r) {
  d <- draw()
  fh <- scfh_sae(y ~ x, "vardir", data = d, W = W, K = 1, mse = "none")
  do.call(rbind, lapply(phis, function(phi) {
    f <- suppressWarnings(scfh_sae(y ~ x, "vardir", data = d, W = W, K = 3, phi = phi, mse = "none"))
    data.frame(r = r, phi = phi, ari = ari(f$cluster$cluster, regime), share = share_correct(f$cluster$cluster),
               rel_rmse = rmse(f, d) / rmse(fh, d), admissible = f$admissible)
  }))
}, mc.cores = cores, mc.set.seed = TRUE)
part1 <- do.call(rbind, part1)
s1 <- aggregate(cbind(ari, share, rel_rmse, admissible) ~ phi, part1, mean)
s1$ari_mcse <- aggregate(ari ~ phi, part1, function(v) sd(v) / sqrt(length(v)))$ari
s1$rel_rmse_mcse <- aggregate(rel_rmse ~ phi, part1, function(v) sd(v) / sqrt(length(v)))$rel_rmse
cat("Part 1: K = 3, R =", R1, "\n")
print(signif(s1, 3), row.names = FALSE)

# Part 2
t0 <- Sys.time()
part2 <- parallel::mclapply(seq_len(R2), function(r) {
  d <- draw()
  f <- suppressWarnings(scfh_sae(y ~ x, "vardir", data = d, W = W, K = 3, phi = 0.5, B = B))
  list(err = f$df_est$est - d$theta, mse = f$df_est$mse, B = f$B_success)
}, mc.cores = cores, mc.set.seed = TRUE)
mc_mse <- rowMeans(sapply(part2, `[[`, "err")^2)
boot_mse <- rowMeans(sapply(part2, `[[`, "mse"))
ratio <- boot_mse / mc_mse
cat("\nPart 2: K = 3, phi = 0.5, R =", R2, ", B =", B, "(", format(Sys.time() - t0), ")\n")
cat("median area ratio boot / MC MSE:", signif(median(ratio), 3),
    "; overall (ratio of means):", signif(mean(boot_mse) / mean(mc_mse), 3),
    "; P10-P90:", signif(quantile(ratio, c(0.1, 0.9)), 3), "\n")
cat("by regime (ratio of means):", signif(tapply(boot_mse, regime, mean) / tapply(mc_mse, regime, mean), 3), "\n")
cat("successful replicates: min", min(sapply(part2, `[[`, "B")), "\n")
