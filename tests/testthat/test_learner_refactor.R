test_that("merf_sae() is unchanged by the shared learner engine", {
  # fixture made with the code before .merf_fit() took a learner
  ref <- readRDS(test_path("fixtures", "merf_reference.rds"))
  # ranger's forests depend on the platform's C++ random number generation, so the fixture
  # is reproducible only on the platform and ranger version that made it
  skip_on_cran()
  skip_if(Sys.info()[["sysname"]] != attr(ref, "sysname") || packageVersion("ranger") != attr(ref, "ranger"),
          "fixture made on another platform or ranger version")
  d <- sim_nonlinear_area(D = 60, seed = 11)
  d$N <- withr::with_seed(1, rpois(60, 300))
  d$y[c(4, 40)] <- NA
  pick <- function(f) list(df_est = f$df_est, coefficients = f$coefficients, sigma2_v = f$sigma2_v,
                           iterations = f$iterations, boot_errors = f$boot_errors)
  for (v in c("replace", "correct")) {
    for (a in c("fhrf", "merf")) {
      f <- merf_sae(y ~ x1 + x2, "vardir", data = d, variant = v, algorithm = a, num.trees = 100, B = 5, seed = 3)
      expect_identical(pick(f), ref[[paste(v, a)]])
    }
  }
  f <- merf_sae(y ~ x1 + x2, "vardir", data = d, variant = "correct", bench_weight = "N", num.trees = 100,
                B = 5, seed = 3)
  expect_identical(pick(f), ref[["bench"]])
})
