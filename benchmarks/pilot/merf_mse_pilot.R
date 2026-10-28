# Pilot check of the bootstrap MSE of merf_sae(algorithm = "fhrf"): model-based
# Monte Carlo with fixed covariates, R = 100 replications, D = 100, B = 50.
# Result (28 Oct 2026): technical-report/DISCREPANCIES.md, check S6.
library(flexsae) # installed package: R CMD INSTALL . --preclean
set.seed(20261028)
D <- 100; R <- 100; B <- 50
x <- data.frame(x1 = runif(D, -2, 2), x2 = rnorm(D)); vardir <- runif(D, 0.2, 0.6)
res <- list()
for (v in c("replace", "correct")) {
  err2 <- boot <- matrix(NA, D, R)
  for (r in 1:R) {
    d <- x; d$vardir <- vardir
    theta <- 2 + sin(1.5 * d$x1) + 0.5 * d$x2 + rnorm(D, sd = 0.5)
    d$y <- theta + rnorm(D, sd = sqrt(vardir))
    f <- merf_sae(y ~ x1 + x2, "vardir", data = d, variant = v, num.trees = 200, B = B, seed = r)
    err2[, r] <- (f$df_est$est - theta)^2; boot[, r] <- f$df_est$mse
  }
  emp <- rowMeans(err2); bm <- rowMeans(boot)
  rb <- (bm - emp) / emp
  cat(sprintf("%-8s mean emp MSE %.4f  mean boot MSE %.4f  ratio %.3f  area RB median %.3f [%.3f, %.3f]\n",
              v, mean(emp), mean(bm), mean(bm) / mean(emp), median(rb), quantile(rb, .1), quantile(rb, .9)))
}
