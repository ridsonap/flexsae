test_that("print and summary describe the fit", {
  dat <- sim_nonlinear_area(D = 40)
  dat$y[3] <- NA
  fit <- pspline_sae(y ~ s(x1, k = 5) + x2, vardir = "vardir", domain = "area", data = dat, B = 5, seed = 1)

  out <- capture.output(print(fit))
  expect_match(out[1], "penalized splines")
  expect_true(any(grepl("40 \\(39 sampled, 1 unsampled\\)", out)))
  expect_true(any(grepl("s\\(x1\\)", out)))
  expect_true(any(grepl("B = 5 \\(5 successful\\)", out)))
  expect_false(any(grepl("flexsae_area", out)))

  s <- summary(fit)
  expect_s3_class(s, "summary.flexsae")
  expect_identical(rownames(s$quantiles), c("est", "rse", "direct"))
  expect_output(print(s), "Distribution over areas")

  none <- pspline_sae(y ~ x1, vardir = "vardir", data = dat, mse = "none")
  expect_output(print(none), "MSE: not estimated")
  expect_identical(rownames(summary(none)$quantiles), c("est", "direct"))
})

test_that("coef, fitted and residuals extract the results", {
  dat <- sim_nonlinear_area(D = 40)
  dat$y[3] <- NA
  fit <- pspline_sae(y ~ s(x1, k = 5) + x2, vardir = "vardir", domain = "area", data = dat, mse = "none")
  df <- fit$df_est

  expect_named(coef(fit), c("(Intercept)", "x2"))
  expect_equal(unname(fitted(fit)), df$est)
  expect_named(fitted(fit), dat$area)
  r <- residuals(fit)
  expect_equal(unname(r), df$direct - df$est)
  expect_true(is.na(r[3]))
  expect_equal(unname(residuals(fit, type = "pearson")),
               (df$direct - df$synthetic) / sqrt(fit$sigma2_v + df$vardir))
  expect_error(residuals(fit, type = "deviance"), "should be one of")
})
