wsum_gap <- function(fit, w, g = rep("All", length(w))) {
  de <- fit$df_est
  s <- de$sampled
  tapply(w[s] * (de$direct[s] - de$est[s]), g[s], sum)
}

test_that("augmented linear model equals sae::eblupFH with the z column (Wang et al., 2008, eq. 23)", {
  skip_if_not_installed("sae")
  milk <- milk_data()
  milk$N <- withr::with_seed(1, rpois(nrow(milk), 300))
  fit <- pspline_sae(yi ~ as.factor(MajorArea), "var", data = milk, bench_weight = "N", mse = "none")
  milk$z <- milk$N * milk$var
  ref <- sae::eblupFH(yi ~ as.factor(MajorArea) + z, var, data = milk, PRECISION = 1e-10)
  expect_equal(fit$sigma2_v, ref$fit$refvar, tolerance = 1e-5)
  expect_equal(fit$df_est$est, unname(drop(ref$eblup)), tolerance = 1e-6)
  expect_true(fit$self_benchmark)
  expect_true("bench(All)" %in% rownames(fit$coefficients))
  expect_lt(abs(wsum_gap(fit, milk$N)), 1e-8)
})

test_that("pspline_sae: restriction holds with smooths, groups and unsampled areas", {
  d <- sim_nonlinear_area(D = 80, seed = 3)
  d$N <- withr::with_seed(1, rpois(80, 500))
  d$g <- rep(c("a", "b"), 40)
  d$y[c(5, 50)] <- NA
  for (method in c("REML", "ML")) {
    fit <- pspline_sae(y ~ s(x1) + x2, "vardir", data = d, method = method,
                       bench_weight = "N", bench_group = "g", mse = "none")
    expect_lt(max(abs(wsum_gap(fit, d$N, d$g))), 1e-7)
    expect_true(all(c("bench(a)", "bench(b)") %in% rownames(fit$coefficients)))
  }
  # unsampled areas get z = 0: synthetic prediction without the augmentation term
  plain <- pspline_sae(y ~ s(x1) + x2, "vardir", data = d, mse = "none")
  expect_false(isTRUE(all.equal(fit$df_est$est, plain$df_est$est)))
  # a weighted-mean target of the direct estimates needs no adjustment afterwards
  bm <- benchmark_flexsae(fit, target = "direct", weight = "N", group = "g", areas = "sampled",
                          method = "difference", mse = "none")
  expect_lt(max(abs(bm$adjustment)), 1e-8)
})

test_that("merf_sae: restriction holds for both variants and algorithms", {
  d <- sim_nonlinear_area(D = 80, seed = 4)
  d$N <- withr::with_seed(1, rpois(80, 500))
  d$g <- rep(c("a", "b"), 40)
  d$y[3] <- NA
  for (variant in c("replace", "correct")) {
    for (algorithm in c("fhrf", "merf")) {
      fit <- merf_sae(y ~ x1 + x2, "vardir", data = d, variant = variant, algorithm = algorithm,
                      bench_weight = "N", bench_group = "g", num.trees = 100, mse = "none", seed = 1)
      expect_lt(max(abs(wsum_gap(fit, d$N, d$g))), 1e-8)
      expect_true(fit$self_benchmark)
      expect_true(all(c("bench(a)", "bench(b)") %in% rownames(fit$coefficients)))
      de <- fit$df_est
      expect_equal(de$est[de$sampled], (de$synthetic + de$gamma * (de$direct - de$synthetic))[de$sampled])
    }
  }
})

test_that("only the relative size of the weights matters", {
  d <- sim_nonlinear_area(D = 50, seed = 5)
  d$N <- withr::with_seed(1, rpois(50, 500))
  a <- pspline_sae(y ~ x1 + x2, "vardir", data = d, bench_weight = "N", mse = "none")
  b <- pspline_sae(y ~ x1 + x2, "vardir", data = d, bench_weight = d$N / sum(d$N), mse = "none")
  expect_equal(a$df_est$est, b$df_est$est, tolerance = 1e-6)
})

test_that("the bootstrap refits the augmented model", {
  d <- sim_nonlinear_area(D = 40, seed = 6)
  d$N <- withr::with_seed(1, rpois(40, 500))
  fit <- pspline_sae(y ~ x1 + x2, "vardir", data = d, bench_weight = "N", B = 10, seed = 2)
  est_b <- fit$boot_errors + fit$boot_theta
  # every replicate satisfies the restriction for its own direct estimates
  gaps <- colSums(d$N * (fit$boot_y - est_b))
  expect_lt(max(abs(gaps)), 1e-6)
  mf <- merf_sae(y ~ x1 + x2, "vardir", data = d, bench_weight = "N", num.trees = 50, B = 5, seed = 2)
  gaps <- colSums(d$N * (mf$boot_y - (mf$boot_errors + mf$boot_theta)))
  expect_lt(max(abs(gaps)), 1e-8)
})

test_that("z in the span of the covariates still gives the restriction", {
  d <- sim_nonlinear_area(D = 30, seed = 7)
  d$vardir <- 0.4 # with equal weights, z is constant like the intercept
  w <- rep(1, 30)
  expect_no_warning(fit <- pspline_sae(y ~ x1 + x2, "vardir", data = d, bench_weight = w, mse = "none"))
  expect_lt(abs(wsum_gap(fit, w)), 1e-8)
  mf <- merf_sae(y ~ x1 + x2, "vardir", data = d, variant = "correct", bench_weight = w,
                 num.trees = 50, mse = "none", seed = 1)
  expect_lt(abs(wsum_gap(mf, w)), 1e-8)
})

test_that(".check_self_benchmark() warns when the restriction fails", {
  sb <- list(weight = rep(1, 4), group = rep("All", 4), groups = "All")
  expect_no_warning(.check_self_benchmark(sb, c(1, 2, 3, 4), c(1.5, 1.5, 3.5, 3.5), rep(TRUE, 4)))
  expect_warning(.check_self_benchmark(sb, c(1, 2, 3, 4), c(2, 2, 3, 4), rep(TRUE, 4)), "does not hold")
})

test_that("invalid self-benchmarking input gives clear errors", {
  d <- sim_nonlinear_area(D = 30, seed = 8)
  expect_error(pspline_sae(y ~ x1, "vardir", data = d, bench_group = rep(1, 30)), "needs")
  expect_error(pspline_sae(y ~ x1, "vardir", data = d, bench_weight = rep(-1, 30)), "positive")
  expect_error(merf_sae(y ~ x1, "vardir", data = d, bench_weight = rep(1, 30),
                        bench_group = c(NA, rep(1, 29))), "missing")
  d$.flexsae_z1 <- 1
  expect_error(pspline_sae(y ~ x1, "vardir", data = d), "reserved")
})
