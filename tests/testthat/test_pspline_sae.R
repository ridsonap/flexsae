# ---- agreement with the linear Fay-Herriot model ------------------------------

test_that("without smooth terms the fit equals sae::eblupFH (REML and ML)", {
  skip_if_not_installed("sae")
  milk <- milk_data()
  for (m in c("REML", "ML")) {
    ref <- sae::eblupFH(yi ~ as.factor(MajorArea), vardir = var, data = milk, method = m)
    fit <- pspline_sae(yi ~ as.factor(MajorArea), vardir = "var", data = milk, method = m, mse = "none")
    expect_equal(fit$sigma2_v, ref$fit$refvar, tolerance = 1e-5)
    expect_equal(unname(coef(fit)), ref$fit$estcoef$beta, tolerance = 1e-6)
    expect_equal(fit$df_est$est, as.numeric(ref$eblup), tolerance = 1e-6)
  }
})

test_that("without smooth terms the fit equals fastsae::eblup_fh", {
  milk <- milk_data()
  milk$dom <- milk$SmallArea
  fit <- pspline_sae(yi ~ as.factor(MajorArea), vardir = "var", domain = "dom", data = milk, mse = "none")

  skip_if_not_installed("fastsae")
  fs <- fastsae::eblup_fh(yi ~ as.factor(MajorArea), vardir = "var", data = milk,
                          precision = 1e-10, print_result = FALSE)
  expect_equal(fit$sigma2_v, fs$random_effect_var, tolerance = 1e-5)
  expect_equal(fit$df_est$est, fs$df_eblup$eblup, tolerance = 1e-6)
})

test_that("without smooth terms the fit equals emdi::fh up to emdi's optimizer tolerance", {
  skip_if_not_installed("emdi")
  milk <- milk_data()
  milk$dom <- milk$SmallArea
  fit <- pspline_sae(yi ~ as.factor(MajorArea), vardir = "var", domain = "dom", data = milk, mse = "none")
  # emdi maximizes the REML likelihood with stats::optimize() at its default
  # tolerance (about 1.2e-4); its `tol` argument does not apply to this case,
  # so sigma2_v differs by about 1.5e-4 relative (DISCREPANCIES.md, check S2).
  e <- suppressMessages(emdi::fh(
    fixed = yi ~ as.factor(MajorArea), vardir = "var", combined_data = milk,
    domains = "dom", method = "reml", MSE = FALSE
  ))
  expect_equal(fit$sigma2_v, e$model$variance, tolerance = 5e-4)
  expect_equal(fit$df_est$est, e$ind$FH, tolerance = 1e-4)
})

test_that("unsampled areas get the synthetic estimate of the fit on the sampled areas", {
  skip_if_not_installed("sae")
  milk <- milk_data()
  out <- c(2, 10, 30)
  milk_u <- milk
  milk_u$yi[out] <- NA
  fit <- pspline_sae(yi ~ as.factor(MajorArea), vardir = "var", data = milk_u, mse = "none")
  ref <- sae::eblupFH(yi ~ as.factor(MajorArea), vardir = var, data = milk[-out, ], method = "REML")
  X <- stats::model.matrix(~ as.factor(MajorArea), milk)

  expect_equal(fit$df_est$sampled, !seq_len(nrow(milk)) %in% out)
  expect_equal(fit$df_est$est[out], unname(drop(X[out, ] %*% ref$fit$estcoef$beta)), tolerance = 1e-6)
  expect_equal(fit$df_est$random_effect[out], rep(0, 3))
  expect_true(all(is.na(fit$df_est$gamma[out])))
  expect_equal(fit$df_est$est[-out], as.numeric(ref$eblup), tolerance = 1e-6)
})

# ---- the non-linear model ---------------------------------------------------------

test_that("random effects are the shrunken residuals of the synthetic part", {
  dat <- sim_nonlinear_area()
  fit <- pspline_sae(y ~ s(x1) + x2, vardir = "vardir", data = dat, mse = "none")
  df <- fit$df_est
  expect_equal(df$random_effect, df$gamma * (df$direct - df$synthetic), tolerance = 1e-8)
  expect_equal(df$est, df$synthetic + df$random_effect)
  expect_equal(df$gamma, fit$sigma2_v / (fit$sigma2_v + dat$vardir))
})

test_that("a smooth term beats the linear model when the effect is non-linear", {
  dat <- sim_nonlinear_area(D = 200, seed = 7)
  lin <- pspline_sae(y ~ x1 + x2, vardir = "vardir", data = dat, mse = "none")
  spl <- pspline_sae(y ~ s(x1) + x2, vardir = "vardir", data = dat, mse = "none")
  rmse <- function(f) sqrt(mean((f$df_est$est - dat$theta)^2))
  expect_lt(rmse(spl), rmse(lin))
  expect_gt(spl$smooth[["edf"]], 2)
})

test_that("the formula may transform the response and use mgcv smooth options", {
  dat <- sim_nonlinear_area()
  dat$ypos <- exp(dat$y / 4)
  fit <- pspline_sae(log(ypos) ~ s(x1, k = 5, bs = "ps") + x2, vardir = "vardir", data = dat, mse = "none")
  expect_equal(fit$df_est$direct, dat$y / 4)
  expect_identical(rownames(fit$smooth), "s(x1)")
})

# ---- bootstrap MSE ------------------------------------------------------------------

test_that("the bootstrap MSE is reproducible and leaves the session seed alone", {
  dat <- sim_nonlinear_area(D = 40)
  set.seed(99)
  before <- .Random.seed
  f1 <- pspline_sae(y ~ s(x1, k = 5) + x2, vardir = "vardir", data = dat, B = 10, seed = 3)
  expect_identical(.Random.seed, before)
  f2 <- pspline_sae(y ~ s(x1, k = 5) + x2, vardir = "vardir", data = dat, B = 10, seed = 3)
  expect_identical(f1$df_est$mse, f2$df_est$mse)
  expect_identical(f1$B_success, 10L)
  expect_true(all(f1$df_est$mse > 0))
  expect_equal(dim(f1$boot_errors), c(40L, 10L))
  expect_equal(f1$df_est$mse, rowMeans(f1$boot_errors^2))
})

test_that("for the linear model the bootstrap MSE agrees with the Prasad-Rao MSE on average", {
  skip_on_cran()
  skip_if_not_installed("sae")
  milk <- milk_data()
  ref <- sae::mseFH(yi ~ as.factor(MajorArea), vardir = var, data = milk, method = "REML")
  fit <- pspline_sae(yi ~ as.factor(MajorArea), vardir = "var", data = milk, B = 400, seed = 1)
  ratio <- mean(fit$df_est$mse) / mean(ref$mse)
  expect_gt(ratio, 0.85)
  expect_lt(ratio, 1.15)
})

test_that("intervals and rse follow from the MSE; mse = 'none' keeps the columns", {
  dat <- sim_nonlinear_area(D = 40)
  fit <- pspline_sae(y ~ s(x1, k = 5) + x2, vardir = "vardir", data = dat, B = 5, seed = 1, conf_level = 0.9)
  df <- fit$df_est
  z <- qnorm(0.95)
  expect_equal(df$ci_upper - df$est, z * sqrt(df$mse))
  expect_equal(df$rse, 100 * sqrt(df$mse) / abs(df$est))

  none <- pspline_sae(y ~ s(x1, k = 5) + x2, vardir = "vardir", data = dat, mse = "none")
  expect_named(none$df_est, names(df))
  expect_true(all(is.na(none$df_est$mse)))
  expect_true(is.na(none$B))
})

# ---- input checks -------------------------------------------------------------------

test_that("invalid input is rejected before fitting", {
  dat <- sim_nonlinear_area(D = 30)
  f <- y ~ s(x1, k = 5) + x2
  expect_error(pspline_sae(f, vardir = "nope", data = dat), "not found")
  expect_error(pspline_sae(y ~ s(z), vardir = "vardir", data = dat), "not found")
  expect_error(pspline_sae(~ x1, vardir = "vardir", data = dat), "two-sided")
  expect_error(pspline_sae(f, vardir = "vardir", data = as.list(dat)), "data frame")

  bad <- dat
  bad$vardir[3] <- 0
  expect_error(pspline_sae(f, vardir = "vardir", data = bad), "positive")
  bad <- dat
  bad$x1[4] <- NA
  expect_error(pspline_sae(f, vardir = "vardir", data = bad), "missing values")
  bad <- dat
  bad$y[-(1:2)] <- NA
  expect_error(pspline_sae(f, vardir = "vardir", data = bad), "At least 3")
  bad <- dat
  bad$.flexsae_area <- 1
  expect_error(pspline_sae(f, vardir = "vardir", data = bad), "reserved")

  expect_error(pspline_sae(f, vardir = "vardir", data = dat, B = 0), "whole number")
  expect_error(pspline_sae(f, vardir = "vardir", data = dat, B = 2.5), "whole number")
  expect_error(pspline_sae(f, vardir = "vardir", data = dat, seed = "a"), "seed")
  expect_error(pspline_sae(f, vardir = "vardir", data = dat, conf_level = 1), "conf_level")
  expect_error(pspline_sae(f, vardir = "vardir", data = dat, method = "GCV"), "should be one of")
  expect_error(pspline_sae(y ~ s(x1, k = 40), vardir = "vardir", data = dat, mse = "none"), "could not be fitted")
})

test_that("vardir and domain accept column names, formulas and vectors", {
  dat <- sim_nonlinear_area(D = 30)
  a <- pspline_sae(y ~ x1, vardir = "vardir", domain = "area", data = dat, mse = "none")
  b <- pspline_sae(y ~ x1, vardir = ~vardir, domain = ~area, data = dat, mse = "none")
  c <- pspline_sae(y ~ x1, vardir = dat$vardir, domain = dat$area, data = dat, mse = "none")
  expect_equal(a$df_est, b$df_est)
  expect_equal(a$df_est, c$df_est)
  expect_identical(a$df_est$domain, dat$area)

  dup <- dat
  dup$area[2] <- dup$area[1]
  expect_warning(pspline_sae(y ~ x1, vardir = "vardir", domain = "area", data = dup, mse = "none"), "duplicated")
})

test_that("unsampled areas need covariate levels seen in the sampled areas", {
  dat <- sim_nonlinear_area(D = 30)
  dat$g <- rep(c("a", "b", "c"), 10)
  dat$y[dat$g == "c"] <- NA
  expect_error(pspline_sae(y ~ x1 + g, vardir = "vardir", data = dat, mse = "none"), "Prediction for the areas failed")
})
