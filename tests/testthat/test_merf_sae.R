test_that(".fit_fh_lm() reproduces sae::eblupFH (REML and ML)", {
  skip_if_not_installed("sae")
  milk <- milk_data()
  X <- model.matrix(~ as.factor(MajorArea), milk)
  for (method in c("REML", "ML")) {
    ours <- .fit_fh_lm(milk$yi, X, milk$var, method)
    ref <- sae::eblupFH(yi ~ as.factor(MajorArea), var, method = method, data = milk)
    expect_equal(ours$sigma2_v, ref$fit$refvar, tolerance = 1e-5)
    expect_equal(unname(ours$beta), ref$fit$estcoef$beta, tolerance = 1e-6)
    expect_equal(ours$coefficients[["Std. Error"]], ref$fit$estcoef$std.error, tolerance = 1e-5)
    expect_true(ours$converged)
  }
})

test_that(".fit_fh_lm() without fixed effects maximizes the likelihood, REML = ML", {
  milk <- milk_data()
  y <- milk$yi - 0.9
  ll <- function(s) -0.5 * sum(log(s + milk$var) + y^2 / (s + milk$var))
  ref <- stats::optimize(ll, c(0, 1), maximum = TRUE, tol = 1e-12)$maximum
  X0 <- matrix(numeric(0), length(y), 0)
  expect_equal(.fit_fh_lm(y, X0, milk$var, "ML")$sigma2_v, ref, tolerance = 1e-6)
  expect_equal(.fit_fh_lm(y, X0, milk$var, "REML")$sigma2_v, ref, tolerance = 1e-6)
  expect_equal(nrow(.fit_fh_lm(y, X0, milk$var, "ML")$coefficients), 0)
})

test_that(".fit_fh_lm() returns 0 at the boundary", {
  withr::local_seed(3)
  vardir <- rep(1, 30)
  y <- rnorm(30, sd = 0.3) # var(y) far below vardir
  expect_equal(.fit_fh_lm(y, matrix(1, 30, 1), vardir, "REML")$sigma2_v, 0)
})

test_that("the EBLUP step with a linear fixed part is the linear Fay-Herriot EBLUP", {
  skip_if_not_installed("sae")
  milk <- milk_data()
  ref <- sae::eblupFH(yi ~ as.factor(MajorArea), var, method = "REML", data = milk)
  X <- model.matrix(~ as.factor(MajorArea), milk)
  fh <- .fit_fh_lm(milk$yi, X, milk$var, "REML")
  out <- .fhrf_combine(drop(X %*% fh$beta), milk$yi, fh$sigma2_v, milk$var, rep(TRUE, nrow(milk)))
  expect_equal(out$est, drop(ref$eblup), tolerance = 1e-6)
})

test_that("fhrf: sigma2_v and coefficients come from the linear Fay-Herriot fit", {
  skip_if_not_installed("sae")
  d <- sim_nonlinear_area(D = 60, seed = 4)
  ref <- sae::eblupFH(y ~ x1 + x2, vardir, method = "REML", data = d)
  fit <- merf_sae(y ~ x1 + x2, "vardir", data = d, variant = "correct",
                  num.trees = 100, mse = "none", seed = 1)
  expect_equal(fit$sigma2_v, ref$fit$refvar, tolerance = 1e-5)
  expect_equal(fit$coefficients$Estimate, ref$fit$estcoef$beta, tolerance = 1e-6)
  expect_equal(fit$iterations, 1L)

  rep <- merf_sae(y ~ x1 + x2, "vardir", data = d, num.trees = 100, mse = "none", seed = 1)
  expect_equal(rep$sigma2_v, ref$fit$refvar, tolerance = 1e-5)
  expect_equal(nrow(rep$coefficients), 0)
  expect_length(coef(rep), 0)
})

test_that("fixed part uses out-of-bag predictions and the predictor follows eqs. (10) and (14)", {
  d <- sim_nonlinear_area(D = 60, seed = 5)
  d$y[c(3, 30)] <- NA
  s <- !is.na(d$y)
  for (variant in c("replace", "correct")) {
    fit <- merf_sae(y ~ x1 + x2, "vardir", data = d, variant = variant,
                    num.trees = 100, mse = "none", seed = 2)
    de <- fit$df_est
    lin <- if (variant == "correct") drop(model.matrix(~ x1 + x2, d) %*% fit$fit$beta) else 0
    f <- unname(de$synthetic - lin)
    expect_equal(f[s], fit$rf$predictions)
    expect_equal(f[!s], predict(fit$rf, data = d[!s, c("x1", "x2")])$predictions)
    expect_equal(de$gamma[s], fit$sigma2_v / (fit$sigma2_v + d$vardir[s]))
    expect_equal(de$est[s], de$synthetic[s] + de$gamma[s] * (d$y[s] - de$synthetic[s]))
    expect_equal(de$est[!s], de$synthetic[!s])
    expect_equal(de$random_effect[!s], c(0, 0))
    expect_true(all(is.na(de$gamma[!s])))
  }
})

test_that("the forest improves on the linear model when the effect is non-linear", {
  d <- sim_nonlinear_area(D = 300, seed = 2)
  rmse <- function(est) sqrt(mean((est - d$theta)^2))
  lin <- pspline_sae(y ~ x1 + x2, "vardir", data = d, mse = "none")
  rf <- merf_sae(y ~ x1 + x2, "vardir", data = d, mse = "none", seed = 1)
  expect_lt(rmse(rf$df_est$est), rmse(lin$df_est$est))
})

test_that("merf algorithm iterates, re-estimates sigma2_v and reports non-convergence", {
  d <- sim_nonlinear_area(D = 150, seed = 2)
  for (variant in c("replace", "correct")) {
    fit <- merf_sae(y ~ x1 + x2, "vardir", data = d, variant = variant,
                    algorithm = "merf", num.trees = 200, mse = "none", seed = 1)
    expect_true(fit$convergence)
    expect_gt(fit$iterations, 1)
    expect_lte(fit$iterations, 25)
    fhrf <- merf_sae(y ~ x1 + x2, "vardir", data = d, variant = variant,
                     num.trees = 200, mse = "none", seed = 1)
    # the linear model puts the missed non-linearity into sigma2_v
    expect_lt(fit$sigma2_v, fhrf$sigma2_v)
  }
  expect_warning(
    fit1 <- merf_sae(y ~ x1 + x2, "vardir", data = d, algorithm = "merf", maxit = 1,
                     num.trees = 50, mse = "none", seed = 1),
    "did not converge in 1 iteration"
  )
  expect_false(fit1$convergence)
  expect_output(print(fit1), "did not fully converge")
})

test_that("seed makes the fit reproducible and restores the session stream", {
  d <- sim_nonlinear_area(D = 50, seed = 6)
  withr::local_seed(99)
  before <- .Random.seed
  a <- merf_sae(y ~ x1 + x2, "vardir", data = d, num.trees = 50, B = 5, seed = 10)
  expect_identical(.Random.seed, before)
  b <- merf_sae(y ~ x1 + x2, "vardir", data = d, num.trees = 50, B = 5, seed = 10)
  expect_identical(a$df_est, b$df_est)
  c <- merf_sae(y ~ x1 + x2, "vardir", data = d, num.trees = 50, B = 5, seed = 11)
  expect_false(identical(a$df_est$est, c$df_est$est))
})

test_that("bootstrap MSE has the shared schema", {
  d <- sim_nonlinear_area(D = 50, seed = 7)
  d$y[2] <- NA
  fit <- merf_sae(y ~ x1 + x2, "vardir", data = d, num.trees = 50, B = 8, seed = 3)
  ps <- pspline_sae(y ~ x1 + x2, "vardir", data = d, mse = "none")
  expect_identical(names(fit$df_est), names(ps$df_est))
  expect_s3_class(fit, c("flexsae_merf", "flexsae"))
  expect_equal(fit$B_success, 8L)
  expect_equal(dim(fit$boot_errors), c(50L, 8L))
  expect_true(all(is.finite(fit$df_est$mse) & fit$df_est$mse > 0))
  expect_equal(fit$df_est$mse, rowMeans(fit$boot_errors^2))

  none <- merf_sae(y ~ x1 + x2, "vardir", data = d, num.trees = 50, mse = "none", seed = 3)
  expect_identical(names(none), names(fit))
  expect_true(all(is.na(none$df_est$mse)))
  expect_identical(none$df_est$est, fit$df_est$est)
})

test_that("factor and character covariates are accepted", {
  d <- sim_nonlinear_area(D = 60, seed = 8)
  d$g <- rep(c("a", "b", "c"), 20)
  fit <- merf_sae(y ~ x1 + g, "vardir", data = d, variant = "correct",
                  num.trees = 50, mse = "none", seed = 1)
  expect_equal(nrow(fit$coefficients), 4) # intercept, x1, gb, gc
  expect_true(all(is.finite(fit$df_est$est)))
})

test_that("invalid input gives clear errors", {
  d <- sim_nonlinear_area(D = 40, seed = 9)
  expect_error(merf_sae(y ~ s(x1), "vardir", data = d), "smooth terms")
  expect_error(merf_sae(y ~ 1, "vardir", data = d), "at least one covariate")
  expect_error(merf_sae(y ~ x1, "vardir", data = d, mtry = 2), "cannot exceed")
  expect_error(merf_sae(y ~ x1, "vardir", data = d, num.trees = 0), "num.trees")
  expect_error(merf_sae(y ~ x1, "vardir", data = d, min.node.size = 1.5), "min.node.size")
  expect_error(merf_sae(y ~ x1, "vardir", data = d, tol = -1), "tol")
  expect_error(merf_sae(y ~ x1, "vardir", data = d, variant = "other"), "should be one of")
  expect_error(merf_sae(y ~ x1, "vardir", data = d, num.trees = 2, mse = "none", seed = 1),
               "never out of bag")
  d$x3 <- 2 * d$x1
  expect_error(merf_sae(y ~ x1 + x3, "vardir", data = d, mse = "none"), "rank deficient")
})
