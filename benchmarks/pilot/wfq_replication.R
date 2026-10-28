# Pilot replication of the simulation of Wang, Fuller and Qu (2008), Section 5,
# with pspline_sae() (linear Fay-Herriot, REML) and self-benchmarking.
# Differences from the paper: sigma2_e / m_i is known (the paper estimates it
# from unit data) and sigma2_b is estimated by REML (the paper uses Wang and
# Fuller, 2003). Result: technical-report/DISCREPANCIES.md, check S8.
library(flexsae) # installed package: R CMD INSTALL . --preclean
pop <- c(33640, 21160, 19360, 16000, 12250, 12250, 11560, 10240, 8410, 8410, 7840, 7290, 6250,
         6250, 5760, 5760, 5760, 5290, 5290, 5290, 4840, 4410, 4410, 4410, 4000, 4000, 3610,
         3240, 3240, 2890, 2890, 2560, 2560, 2250, 1960, 1690, 1690, 1690, 1210, 1210, 1210,
         1210, 1000, 810, 810, 810, 640, 640, 640, 490) / 1000 # millions (Table 1)
m <- c(58, 46, 44, 40, 35, 35, 34, 32, 29, 29, 28, 27, 25, 25, 24, 24, 24, 23, 23, 23, 22, 21,
       21, 21, 20, 20, 19, 18, 18, 17, 17, 16, 16, 15, 14, 13, 13, 13, 11, 11, 11, 11, 10, 9, 9,
       9, 8, 8, 8, 7)
z <- pop^0.2 - mean(pop^0.2)
omega <- pop / sum(pop)
psi <- 16 / m
gam_known <- 1 / (1 + psi) # sigma2_b = 1
phi <- 1 / (gam_known * psi) # eq. (34)
set.seed(20081)
R <- 2000
out <- matrix(NA, R, 6, dimnames = list(NULL, c("Q_A", "Q_AUG_A", "Q_B", "Q_AUG_B", "bias_A", "bias_AUG_A")))
err_A <- err_AUGA <- matrix(NA, R, 50)
for (r in seq_len(R)) {
  y <- 6 + 3 * z + rnorm(50) # theta_i (model B)
  d <- data.frame(Y = y + rnorm(50, sd = sqrt(psi)), z = z, psi = psi, omega = omega)
  fA <- pspline_sae(Y ~ 1, "psi", data = d, mse = "none")$df_est$est
  fAa <- pspline_sae(Y ~ 1, "psi", data = d, bench_weight = "omega", mse = "none")$df_est$est
  fB <- pspline_sae(Y ~ z, "psi", data = d, mse = "none")$df_est$est
  fBa <- pspline_sae(Y ~ z, "psi", data = d, bench_weight = "omega", mse = "none")$df_est$est
  Q <- function(e) 0.02 * sum(phi * (e - y)^2)
  out[r, ] <- c(Q(fA), Q(fAa), Q(fB), Q(fBa), sum(omega * (d$Y - fA)), sum(omega * (d$Y - fAa)))
  err_A[r, ] <- fA - y
  err_AUGA[r, ] <- fAa - y
}
res <- colMeans(out)
mcse <- apply(out, 2, sd) / sqrt(R)
print(round(rbind(mean = res, mcse = mcse), 4))
norm_bias <- function(e) colMeans(e) / sqrt(gam_known * psi)
states <- c(1, 2, 14, 26, 38, 50)
print(round(rbind(EBLUP_A = norm_bias(err_A)[states], AUG2_A = norm_bias(err_AUGA)[states]), 3))
