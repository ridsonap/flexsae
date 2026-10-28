bm_vector_data <- function() {
  withr::with_seed(3, list(
    est = runif(30, 0.2, 0.8),
    w = rpois(30, 200),
    g = rep(c("a", "b", "c"), 10),
    m = runif(30, 0.001, 0.01)
  ))
}

test_that("rules agree with fastsae::benchmark_sae (E4)", {
  skip_if_not_installed("fastsae", "1.0.0")
  d <- bm_vector_data()
  cmp <- function(...) {
    ours <- benchmark_flexsae(d$est, ...)$benchmarked
    ref <- suppressMessages(fastsae::benchmark_sae(d$est, ...))$benchmarked
    max(abs(ours - ref))
  }
  tg_mean <- c(a = 0.5, b = 0.45, c = 0.55)
  for (method in c("ratio", "difference", "logit")) {
    tol <- if (method == "logit") 1e-8 else 1e-10
    for (type in c("mean", "total")) {
      tg <- if (type == "mean") tg_mean else tg_mean * tapply(d$w, d$g, sum)[names(tg_mean)]
      one <- if (type == "mean") 0.52 else 0.52 * sum(d$w)
      expect_lt(cmp(target = one, weight = d$w, method = method, type = type), tol)
      expect_lt(cmp(target = tg, weight = d$w, group = d$g, method = method, type = type), tol)
    }
    expect_lt(cmp(target = tg_mean, weight = d$w, group = d$g, method = method, national_target = 0.52), tol)
    expect_lt(cmp(weight = d$w, group = d$g, method = method, national_target = 0.52), tol)
  }
})

test_that("optimal agrees with fastsae and saebenchmarking", {
  skip_if_not_installed("fastsae", "1.0.0")
  skip_if_not_installed("saebenchmarking")
  d <- bm_vector_data()
  tg <- c(a = 0.5, b = 0.45, c = 0.55)
  ours <- benchmark_flexsae(d$est, target = tg, weight = d$w, group = d$g, method = "optimal", area_mse = d$m)
  ref <- fastsae:::.benchmark_worker(seq_along(d$est), d$est, target = tg, weight = d$w,
                                     method = "optimal", group = d$g, mse = d$m)
  expect_equal(ours$benchmarked, ref$benchmarked, tolerance = 1e-12)

  direct <- d$est + 0.02
  sb <- saebenchmarking::sae_benchmarking("optimum", direct = direct, weight = d$w,
                                          estimate = d$est, mse = d$m)
  ours <- benchmark_flexsae(d$est, target = sum(d$w * direct) / sum(d$w), weight = d$w,
                            method = "optimal", area_mse = d$m)
  expect_equal(ours$benchmarked, sb$estimates$Estimate, tolerance = 1e-12)
})

test_that("ratio and difference agree with emdi::benchmark (E5)", {
  skip_if_not_installed("emdi")
  milk <- milk_data()
  share <- withr::with_seed(2, runif(nrow(milk)))
  share <- share / sum(share)
  efh <- emdi::fh(yi ~ as.factor(MajorArea), vardir = "var", combined_data = milk,
                  domains = "SmallArea", method = "reml", MSE = FALSE)
  for (pair in list(c("raking", "difference"), c("ratio", "ratio"))) {
    ref <- emdi::benchmark(efh, benchmark = 0.95, share = share, type = pair[1])
    ours <- benchmark_flexsae(efh$ind$FH, target = 0.95, weight = share, method = pair[2])
    expect_equal(ours$benchmarked, ref$FH_Bench, tolerance = 1e-12)
  }
})

test_that("benchmarked aggregates equal the targets for every rule", {
  d <- sim_nonlinear_area(D = 60, seed = 11)
  d$theta <- NULL
  d$y <- plogis(d$y - 2) # proportions, for the logit rule
  d$vardir <- d$vardir / 50
  d$N <- withr::with_seed(1, rpois(60, 300))
  d$g <- rep(c("a", "b", "c"), 20)
  fit <- pspline_sae(y ~ x1 + x2, "vardir", data = d, mse = "none")
  fit$df_est$mse <- withr::with_seed(2, runif(60, 0.001, 0.01)) # for "optimal"
  tg <- c(a = 0.6, b = 0.55, c = 0.62)
  for (method in c("ratio", "difference", "optimal", "logit")) {
    tol <- if (method == "logit") 1e-8 else 1e-10
    for (type in c("mean", "total")) {
      if (method == "logit" && type == "total") next
      t_g <- if (type == "mean") tg else tg * tapply(d$N, d$g, sum)[names(tg)]
      bm <- benchmark_flexsae(fit, target = t_g, weight = "N", group = "g", method = method,
                              type = type, mse = "none")
      v <- attr(bm, "verification")
      expect_lt(max(abs(v$benchmarked - t_g[v$group])), tol)
      a <- if (type == "mean") d$N / ave(d$N, d$g, FUN = sum) else d$N
      expect_lt(max(abs(tapply(a * bm$benchmarked, d$g, sum)[names(tg)] - t_g)), tol)
    }
    bm2 <- benchmark_flexsae(fit, target = tg, weight = "N", group = "g", method = method,
                             national_target = 0.6, mse = "none")
    expect_lt(abs(sum(d$N * bm2$benchmarked) / sum(d$N) - 0.6), tol)
    expect_lt(max(abs(attr(bm2, "verification")$discrepancy)), tol)
    expect_true(attr(bm2, "hierarchical"))
  }
})

test_that("areas = 'sampled' leaves unsampled areas unchanged; target = 'direct'", {
  d <- sim_nonlinear_area(D = 50, seed = 12)
  d$y[c(4, 9)] <- NA
  d$N <- withr::with_seed(1, rpois(50, 300))
  fit <- pspline_sae(y ~ x1 + x2, "vardir", data = d, mse = "none")
  expect_error(benchmark_flexsae(fit, target = "direct", weight = "N", mse = "none"),
               "every area in the aggregate")
  bm <- benchmark_flexsae(fit, target = "direct", weight = "N", areas = "sampled", mse = "none")
  s <- !is.na(d$y)
  expect_identical(bm$benchmarked[!s], bm$original[!s])
  expect_true(all(is.na(bm$target[!s])))
  expect_equal(sum(d$N[s] * bm$benchmarked[s]) / sum(d$N[s]), sum(d$N[s] * d$y[s]) / sum(d$N[s]),
               tolerance = 1e-12)
  expect_identical(attr(bm, "target_source"), "direct")
})

test_that("bootstrap MSE re-benchmarks every replicate", {
  d <- sim_nonlinear_area(D = 40, seed = 13)
  d$N <- withr::with_seed(1, rpois(40, 300))
  fit <- pspline_sae(y ~ x1 + x2, "vardir", data = d, B = 30, seed = 5)
  est_b <- fit$boot_errors + fit$boot_theta
  a <- d$N / sum(d$N)
  for (src in c("known", "direct")) {
    bm <- benchmark_flexsae(fit, target = if (src == "known") 3 else "direct", weight = "N",
                            method = "difference")
    tg <- if (src == "known") colSums(a * fit$boot_theta) else colSums(a * fit$boot_y)
    bench_b <- sweep(est_b, 2, tg - colSums(a * est_b), "+")
    expect_equal(bm$mse_benchmarked, rowMeans((bench_b - fit$boot_theta)^2), tolerance = 1e-12)
    expect_equal(bm$mse, fit$df_est$mse, tolerance = 1e-12)
    expect_equal(bm$mse_ratio, bm$mse_benchmarked / bm$mse)
    expect_equal(attr(bm, "B"), 30L)
    expect_equal(attr(bm, "n_failed"), 0L)
  }

  # ratio: rule applied with the replicate's own aggregate
  bm <- benchmark_flexsae(fit, target = 3, weight = "N", method = "ratio")
  tg <- colSums(a * fit$boot_theta)
  bench_b <- sweep(est_b, 2, tg / colSums(a * est_b), "*")
  expect_equal(bm$mse_benchmarked, rowMeans((bench_b - fit$boot_theta)^2), tolerance = 1e-12)
})

test_that("an estimated target adds its variance and is reproducible", {
  d <- sim_nonlinear_area(D = 40, seed = 14)
  fit <- pspline_sae(y ~ x1 + x2, "vardir", data = d, B = 30, seed = 5)
  known <- benchmark_flexsae(fit, target = 3, method = "difference")
  withr::local_seed(42)
  before <- .Random.seed
  est1 <- benchmark_flexsae(fit, target = 3, method = "difference", target_var = 0.05, seed = 1)
  expect_identical(.Random.seed, before)
  est2 <- benchmark_flexsae(fit, target = 3, method = "difference", target_var = 0.05, seed = 1)
  expect_identical(est1$mse_benchmarked, est2$mse_benchmarked)
  expect_identical(attr(est1, "target_source"), "estimated")
  expect_identical(est1$benchmarked, known$benchmarked)
  # difference rule with one group: every area gets the same target error
  eps <- withr::with_seed(1, rnorm(30, 0, sqrt(0.05)))
  bench_b <- fit$boot_errors + fit$boot_theta
  bench_b <- sweep(bench_b, 2, colMeans(fit$boot_theta) + eps - colMeans(bench_b), "+")
  expect_equal(est1$mse_benchmarked, rowMeans((bench_b - fit$boot_theta)^2), tolerance = 1e-12)
  expect_gt(mean(est1$mse_benchmarked), mean(known$mse_benchmarked))
})

test_that("difference MSE increase matches Steorts and Ghosh (2013) via saebenchmarking (E6)", {
  skip_on_cran()
  skip_if_not_installed("saebenchmarking")
  milk <- milk_data()
  milk$N <- withr::with_seed(1, rpois(nrow(milk), 300))
  fit <- pspline_sae(yi ~ as.factor(MajorArea), "var", data = milk, B = 400, seed = 7)
  bm <- benchmark_flexsae(fit, target = "direct", weight = "N", method = "difference")
  sb <- saebenchmarking::sae_benchmarking("difference", direct = milk$yi, weight = milk$N,
                                          estimate = fit$df_est$est, vardir = milk$var)
  z <- model.matrix(~ as.factor(MajorArea), milk)
  ms <- saebenchmarking::mse_benchmarking(sb, z = z, sigma2_v = fit$sigma2_v)
  pr <- sae::mseFH(yi ~ as.factor(MajorArea), var, data = milk)$mse
  g4 <- mean(ms$mse - pr)
  expect_equal(mean(bm$mse_benchmarked - bm$mse), g4, tolerance = 0.3)
})

test_that("merf_sae fits can be benchmarked", {
  d <- sim_nonlinear_area(D = 40, seed = 15)
  fit <- merf_sae(y ~ x1 + x2, "vardir", data = d, num.trees = 50, B = 10, seed = 1)
  bm <- benchmark_flexsae(fit, target = 2.5)
  expect_equal(mean(bm$benchmarked), 2.5, tolerance = 1e-12)
  expect_true(all(is.finite(bm$mse_benchmarked)))
})

test_that("bounds are flagged, print and summary work", {
  d <- sim_nonlinear_area(D = 30, seed = 16)
  fit <- pspline_sae(y ~ x1 + x2, "vardir", data = d, B = 10, seed = 1)
  expect_warning(bm <- benchmark_flexsae(fit, target = 2.5, bounds = c(0, 3)), "outside")
  expect_identical(bm$out_of_bounds, bm$benchmarked < 0 | bm$benchmarked > 3)
  expect_output(print(bm), "Aggregates by group")
  expect_output(summary(bm), "MSE higher after benchmarking")
  expect_s3_class(bm[1:3, ], "data.frame")
  expect_output(print(bm[1:3, ]), "domain")
  expect_identical(names(bm), c("domain", "group", "weight", "sampled", "original", "benchmarked",
                                "adjustment", "rel_adjustment_pct", "target", "mse",
                                "mse_benchmarked", "mse_ratio", "out_of_bounds"))
})

test_that("invalid input gives clear errors", {
  d <- sim_nonlinear_area(D = 30, seed = 17)
  d$g <- rep(c("a", "b"), 15)
  fit <- pspline_sae(y ~ x1 + x2, "vardir", data = d, mse = "none")
  expect_error(benchmark_flexsae(fit, target = 2), "no bootstrap replicates")
  expect_error(benchmark_flexsae(fit, target = 2, method = "optimal", mse = "none"), "needs the MSE")
  expect_error(benchmark_flexsae(fit, target = 2, weight = rep(-1, 30), mse = "none"), "positive")
  expect_error(benchmark_flexsae(fit, target = c(a = 2), group = "g", mse = "none"), "no value for group")
  expect_error(benchmark_flexsae(fit, mse = "none"), "Give")
  expect_error(benchmark_flexsae(fit, target = "direct", target_var = 1, mse = "none"), "target_var")
  expect_error(benchmark_flexsae(fit, target = 2, foo = 1, mse = "none"), "Unknown argument")
  expect_error(benchmark_flexsae(fit, target = 2, bounds = c(3, 1), mse = "none"), "lower < upper")
  expect_error(benchmark_flexsae(fit, target = 2, method = "logit", mse = "none"), "in \\(0, 1\\)")
  expect_error(benchmark_flexsae(1:3, target = "direct"), "needs a")
  expect_error(benchmark_flexsae(c(1, 2), target = 1, area_mse = 1), "length 2")
})

test_that("two-stage bootstrap MSE: known targets make stage 1 a no-op in the replicates", {
  d <- sim_nonlinear_area(D = 40, seed = 18)
  d$g <- rep(c("a", "b"), 20)
  d$N <- withr::with_seed(1, rpois(40, 300))
  fit <- pspline_sae(y ~ x1 + x2, "vardir", data = d, B = 20, seed = 2)
  tg <- c(a = 2.4, b = 2.2)
  two <- benchmark_flexsae(fit, target = tg, weight = "N", group = "g", method = "difference",
                           national_target = 2.5)
  one <- benchmark_flexsae(fit, target = tg, weight = "N", group = "g", method = "difference")
  expect_equal(two$mse_benchmarked, one$mse_benchmarked, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(two$benchmarked, one$benchmarked)))
  free <- benchmark_flexsae(fit, weight = "N", group = "g", method = "ratio", national_target = 2.5)
  expect_true(all(is.finite(free$mse_benchmarked)))
  expect_output(print(free), "National")
})

test_that("summary() does not print an empty area table", {
  d <- sim_nonlinear_area(D = 20, seed = 19)
  fit <- pspline_sae(y ~ x1 + x2, "vardir", data = d, B = 5, seed = 1)
  out <- utils::capture.output(summary(benchmark_flexsae(fit, target = 2.5)))
  expect_false(any(grepl("0 rows|more areas", out)))
  expect_true(any(grepl("Adjustment", out)))
})
