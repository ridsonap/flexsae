skip_if_not_installed("xgboost")

test_that("cross-fitted predictions do not use the area's own target", {
  d <- sim_nonlinear_area(D = 60, seed = 21)
  xmat <- .boost_matrix(y ~ x1 + x2, d)
  sampled <- rep(TRUE, 60)
  sampled[c(5, 6)] <- FALSE
  folds <- rep_len(1:5, 58)
  params <- list(objective = "reg:squarederror", eta = 0.1, max_depth = 2, subsample = 1,
                 min_child_weight = 1, nthread = 1)
  a <- .fit_boost(xmat, d$y, sampled, params, 30, folds, seed = 1)
  target <- d$y
  i <- which(sampled)[1] # fold 1
  target[i] <- target[i] + 50
  b <- .fit_boost(xmat, target, sampled, params, 30, folds, seed = 1)
  s <- which(sampled)
  expect_identical(b$f[s[folds == 1]], a$f[s[folds == 1]])
  expect_false(isTRUE(all.equal(b$f[s[folds == 2]], a$f[s[folds == 2]])))
  # unsampled areas use the booster of all sampled areas
  expect_false(isTRUE(all.equal(b$f[!sampled], a$f[!sampled])))
  expect_equal(a$f[!sampled], unname(stats::predict(a$fit, xmat[!sampled, ])))
})

test_that("the one-pass fit uses the linear Fay-Herriot variance and the EBLUP form", {
  skip_if_not_installed("sae")
  d <- sim_nonlinear_area(D = 80, seed = 22)
  for (method in c("REML", "ML")) {
    fit <- boost_sae(y ~ x1 + x2, "vardir", data = d, method = method, mse = "none", seed = 1)
    ref <- sae::eblupFH(y ~ x1 + x2, vardir, data = d, method = method, PRECISION = 1e-10)
    expect_equal(fit$sigma2_v, ref$fit$refvar, tolerance = 1e-6)
    de <- fit$df_est
    expect_equal(de$est, de$synthetic + de$gamma * (de$direct - de$synthetic))
    expect_equal(de$gamma, fit$sigma2_v / (fit$sigma2_v + de$vardir))
  }
  fit <- boost_sae(y ~ x1 + x2, "vardir", data = d, variant = "correct", mse = "none", seed = 1)
  ref <- sae::eblupFH(y ~ x1 + x2, vardir, data = d, PRECISION = 1e-10)
  expect_equal(unname(fit$coefficients$Estimate), unname(ref$fit$estcoef$beta), tolerance = 1e-5)
})

test_that("boosting improves on the linear model for a non-linear mean", {
  d <- sim_nonlinear_area(D = 200, seed = 23)
  lin <- pspline_sae(y ~ x1 + x2, "vardir", data = d, mse = "none")
  err <- function(f) mean((f$df_est$est - d$theta)^2)
  for (variant in c("replace", "correct")) {
    fit <- boost_sae(y ~ x1 + x2, "vardir", data = d, variant = variant, mse = "none", seed = 1)
    expect_lt(err(fit), err(lin))
  }
})

test_that("results are reproducible, the stream is restored and the bootstrap keeps nrounds", {
  d <- sim_nonlinear_area(D = 50, seed = 24)
  d$y[c(2, 9)] <- NA
  withr::local_seed(10)
  before <- runif(1)
  withr::local_seed(10)
  a <- boost_sae(y ~ x1 + x2, "vardir", data = d, B = 5, seed = 7)
  expect_identical(runif(1), before)
  b <- boost_sae(y ~ x1 + x2, "vardir", data = d, B = 5, seed = 7)
  expect_identical(a$df_est, b$df_est)
  expect_true(a$boost_args$nrounds_chosen)
  expect_identical(a$B_success, 5L)
  expect_equal(a$df_est$mse, rowMeans(a$boot_errors^2))
  expect_equal(a$df_est$est[!a$df_est$sampled], a$df_est$synthetic[!a$df_est$sampled])
  fixed <- boost_sae(y ~ x1 + x2, "vardir", data = d, nrounds = 15, mse = "none", seed = 7)
  expect_identical(fixed$boost_args$nrounds, 15L)
  expect_false(fixed$boost_args$nrounds_chosen)
  bm <- benchmark_flexsae(a, target = 3)
  expect_equal(mean(bm$benchmarked), 3)
})

test_that("self-benchmarking holds for both variants and algorithms", {
  d <- sim_nonlinear_area(D = 60, seed = 25)
  d$N <- withr::with_seed(1, rpois(60, 400))
  d$g <- rep(c("a", "b"), 30)
  for (variant in c("replace", "correct")) {
    for (algorithm in c("onepass", "iterative")) {
      fit <- boost_sae(y ~ x1 + x2, "vardir", data = d, variant = variant, algorithm = algorithm,
                       bench_weight = "N", bench_group = "g", nrounds = 30, mse = "none", seed = 1)
      de <- fit$df_est
      gap <- tapply(d$N * (de$direct - de$est), d$g, sum)
      expect_lt(max(abs(gap)), 1e-8)
      expect_true(all(c("bench(a)", "bench(b)") %in% rownames(fit$coefficients)))
    }
  }
})

test_that("factors, importance, partial dependence and plots", {
  d <- sim_nonlinear_area(D = 80, seed = 26)
  d$f <- rep(c("u", "v", "w", "z"), 20)
  d$y <- d$y + ifelse(d$f == "w", 1, 0)
  fit <- boost_sae(y ~ x1 + x2 + f, "vardir", data = d, nrounds = 40, B = 5, seed = 1)
  expect_identical(fit$boost_args$columns, c("x1", "x2", "fu", "fv", "fw", "fz"))
  imp <- importance_flexsae(fit)
  expect_setequal(imp$variable, c("x1", "x2", "f"))
  expect_equal(sum(imp$importance), 1)
  expect_identical(attr(imp, "type"), "gain")
  expect_error(importance_flexsae(fit, "permutation"), "gain")
  pd <- partial_flexsae(fit, n = 4, ice = TRUE)
  expect_setequal(pd$level[pd$variable == "f"], c("u", "v", "w", "z"))
  ice <- attr(pd, "ice")
  num <- pd[pd$variable == "x1", ]
  means <- tapply(ice$fit[ice$variable == "x1"], ice$x[ice$variable == "x1"], mean)
  expect_equal(as.numeric(means), num$fit, tolerance = 1e-6)
  for (type in c("importance", "partial", "estimates", "residuals", "rse")) {
    expect_silent(ggplot2::ggplot_build(autoplot(fit, type = type)))
  }
  expect_output(print(fit), "40 boosting rounds")
  expect_error(select_flexsae(fit), "not implemented")
  mf <- merf_sae(y ~ x1 + x2, "vardir", data = d, num.trees = 50, mse = "none", seed = 1)
  expect_error(importance_flexsae(mf, "gain"), "boost_sae")
})

test_that("invalid input gives clear errors", {
  d <- sim_nonlinear_area(D = 40, seed = 27)
  expect_error(boost_sae(y ~ x1, "vardir", data = d, eta = 0), "eta")
  expect_error(boost_sae(y ~ x1, "vardir", data = d, subsample = 1.5), "subsample")
  expect_error(boost_sae(y ~ x1, "vardir", data = d, max_depth = 0), "max_depth")
  expect_error(boost_sae(y ~ x1, "vardir", data = d, nfolds = 1), "nfolds")
  expect_error(boost_sae(y ~ x1, "vardir", data = d, nfolds = 50), "nfolds")
  expect_error(boost_sae(y ~ x1, "vardir", data = d, nrounds = 2.5), "nrounds")
  expect_error(boost_sae(y ~ 1, "vardir", data = d), "at least one covariate")
  expect_error(boost_sae(y ~ s(x1), "vardir", data = d), "smooth")
})
