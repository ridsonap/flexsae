sel_data <- function(D = 80, seed = 2) {
  withr::with_seed(seed, {
    d <- data.frame(x1 = runif(D, -2, 2), x2 = rnorm(D), x3 = rnorm(D), x4 = rnorm(D))
    d$vardir <- runif(D, 0.2, 0.6)
    d$y <- 2 + 0.8 * d$x1 + 0.5 * d$x2 + rnorm(D, sd = 0.5) + rnorm(D, sd = sqrt(d$vardir))
    d$y2 <- 2 + sin(1.5 * d$x1) + 0.5 * d$x2 + rnorm(D, sd = 0.5) + rnorm(D, sd = sqrt(d$vardir))
    d
  })
}

test_that("linear stepwise selection equals fastsae::step_sae", {
  skip_if_not_installed("fastsae", "1.0.0")
  d <- sel_data()
  fit <- pspline_sae(y ~ x1 + x2 + x3 + x4, "vardir", data = d, mse = "none")
  ref_fit <- fastsae::eblup_fh(y ~ x1 + x2 + x3 + x4, vardir = "vardir", data = d)
  for (criterion in c("BIC", "AIC", "KIC")) {
    ours <- select_flexsae(fit, criterion = criterion, trace = FALSE)
    utils::capture.output(ref <- suppressMessages(
      fastsae::step_sae(ref_fit, criterion = criterion, direction = "backward", trace = FALSE)
    ))
    expect_equal(ours$path[[3]], ref$path[[2]], tolerance = 1e-6)
    expect_setequal(ours$selected, c("x1", "x2"))
  }
  # the marginal ML log-likelihood is the one of fastsae
  ml <- pspline_sae(y ~ x1 + x2 + x3 + x4, "vardir", data = d, method = "ML", mse = "none")
  ref_ml <- fastsae::eblup_fh(y ~ x1 + x2 + x3 + x4, vardir = "vardir", data = d, method = "ML")
  expect_equal(-as.numeric(ml$fit$gcv.ubre), unname(ref_ml$goodness[["loglikelihood"]]), tolerance = 1e-6)
})

test_that("stepwise selection with smooth terms drops noise terms and refits with the object's settings", {
  d <- sel_data()
  fit <- pspline_sae(y2 ~ s(x1) + s(x2) + s(x3) + x4, "vardir", data = d, B = 5, seed = 1)
  sel <- select_flexsae(fit, trace = FALSE)
  expect_setequal(sel$selected, c("s(x1)", "s(x2)"))
  expect_s3_class(sel$model, "flexsae_pspline")
  expect_identical(sel$model$method, "REML")
  expect_identical(sel$model$B_success, 5L)
  expect_output(print(sel), "Selected model")
  # forward from a small model with a scope
  small <- pspline_sae(y2 ~ x4, "vardir", data = d, mse = "none")
  fw <- select_flexsae(small, direction = "forward", scope = ~ s(x1) + s(x2) + x4, trace = FALSE)
  expect_true(all(c("s(x1)", "s(x2)") %in% fw$selected))
  expect_error(select_flexsae(fit, scope = ~ s(x1)), "must contain the terms")
})

test_that("Boruta keeps the relevant covariates of a merf_sae model", {
  skip_if_not_installed("Boruta")
  d <- sel_data()
  fit <- merf_sae(y2 ~ x1 + x2 + x3 + x4, "vardir", data = d, num.trees = 200, mse = "none", seed = 1)
  sel <- select_flexsae(fit, seed = 1)
  expect_identical(sel$method, "boruta")
  expect_true(all(c("x1", "x2") %in% sel$selected))
  expect_s3_class(sel$model, "flexsae_merf")
  expect_identical(sel$model$rf$num.trees, 200)
  expect_error(select_flexsae(fit, method = "step"), "needs a likelihood")
  expect_error(select_flexsae(pspline_sae(y ~ x1, "vardir", data = d, mse = "none"), method = "boruta"),
               "merf_sae")
})

test_that("importance refits the same forest", {
  d <- sel_data()
  fit <- merf_sae(y2 ~ x1 + x2 + x3 + x4, "vardir", data = d, num.trees = 200, mse = "none", seed = 1)
  imp <- importance_flexsae(fit)
  expect_identical(imp$variable[1], "x1")
  expect_identical(attr(imp, "type"), "permutation")
  same <- ranger::ranger(x = d[, c("x1", "x2", "x3", "x4")], y = fit$rf_target, num.trees = 200,
                         seed = fit$rf_seed, num.threads = 1)
  expect_identical(same$predictions, fit$rf$predictions)
  expect_identical(attr(importance_flexsae(fit, "impurity"), "type"), "impurity")
  expect_error(importance_flexsae(pspline_sae(y ~ x1, "vardir", data = d, mse = "none")), "merf_sae")
})

test_that("partial dependence of an additive spline model follows its terms", {
  d <- sel_data()
  fit <- pspline_sae(y2 ~ s(x1) + x2, "vardir", data = d, mse = "none")
  pd <- partial_flexsae(fit, n = 7)
  lin <- pd[pd$variable == "x2", ]
  expect_equal(diff(lin$fit) / diff(lin$x), rep(fit$coefficients["x2", "Estimate"], 6), tolerance = 1e-8)
  sm <- pd[pd$variable == "x1", ]
  term <- mgcv::predict.gam(fit$fit, newdata = data.frame(x1 = sm$x, x2 = 0, .flexsae_area = factor(1)),
                            type = "terms", terms = "s(x1)")
  expect_equal(diff(sm$fit), diff(unname(term[, 1])), tolerance = 1e-8)
  expect_true(all(sm$lower < sm$fit & sm$fit < sm$upper))

  ice <- attr(partial_flexsae(fit, "x1", n = 4, ice = TRUE), "ice")
  expect_equal(nrow(ice), 4 * nrow(d))
  expect_error(partial_flexsae(fit, "zz"), "not covariate")
})

test_that("partial dependence of merf_sae is the mean of the ICE curves", {
  d <- sel_data()
  d$g <- rep(c("a", "b"), 40)
  d$N <- withr::with_seed(1, rpois(80, 300))
  fit <- merf_sae(y2 ~ x1 + g, "vardir", data = d, variant = "correct", bench_weight = "N",
                  num.trees = 100, mse = "none", seed = 1)
  pd <- partial_flexsae(fit, n = 5, ice = TRUE)
  ice <- attr(pd, "ice")
  num <- pd[pd$variable == "x1", ]
  means <- tapply(ice$fit[ice$variable == "x1"], ice$x[ice$variable == "x1"], mean)
  expect_equal(as.numeric(means), num$fit, tolerance = 1e-12)
  expect_true(all(is.na(pd$se)))
  expect_setequal(pd$level[pd$variable == "g"], c("a", "b"))
  # at the observed covariates the fixed part with ordinary forest predictions
  each <- .fixed_part_fun(fit)(cbind(d, .flexsae_z1 = fit$bench$Z[, 1]))$each
  expect_equal(length(each), nrow(d))
})

test_that("autoplot draws every type", {
  d <- sel_data()
  d$N <- withr::with_seed(1, rpois(80, 300))
  d$y2[c(3, 7)] <- NA
  fit <- pspline_sae(y2 ~ s(x1) + x2, "vardir", data = d, B = 10, seed = 1)
  for (type in c("estimates", "scatter", "rse", "residuals", "qq", "caterpillar", "partial")) {
    p <- autoplot(fit, type = type)
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
  mf <- merf_sae(y2 ~ x1 + x2, "vardir", data = d, num.trees = 50, mse = "none", seed = 1)
  expect_s3_class(autoplot(mf, type = "importance"), "ggplot")
  expect_s3_class(autoplot(mf, type = "partial", ice = TRUE, n = 5), "ggplot")
  expect_error(autoplot(mf, type = "rse"), "no MSE")
  bm <- benchmark_flexsae(fit, target = 2, weight = "N")
  expect_s3_class(autoplot(bm), "ggplot")
  expect_s3_class(autoplot(bm, type = "mse"), "ggplot")
  expect_error(autoplot(benchmark_flexsae(mf, target = 2, mse = "none"), type = "mse"), "no MSE")
})
