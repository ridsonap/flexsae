test_that("K = 1 is the Fay-Herriot model of sae and emdi", {
  skip_if_not_installed("sae")
  milk <- milk_data()
  W <- diag(0, nrow(milk))
  for (method in c("REML", "ML")) {
    fit <- scfh_sae(yi ~ as.factor(MajorArea), "var", data = milk, W = W, K = 1, method = method, mse = "none")
    ref <- sae::eblupFH(yi ~ as.factor(MajorArea), var, data = milk, method = method, PRECISION = 1e-10)
    expect_equal(fit$sigma2_v, ref$fit$refvar, tolerance = 1e-6)
    expect_equal(fit$df_est$est, unname(drop(ref$eblup)), tolerance = 1e-6)
    expect_equal(unname(fit$coefficients$Estimate), unname(ref$fit$estcoef$beta), tolerance = 1e-6)
    if (method == "ML") expect_equal(fit$loglik, unname(ref$fit$goodness[["loglike"]]), tolerance = 1e-6)
  }
  skip_if_not_installed("emdi")
  milk$dom <- seq_len(nrow(milk))
  fit <- scfh_sae(yi ~ as.factor(MajorArea), "var", data = milk, W = W, K = 1, mse = "none")
  ref <- suppressMessages(emdi::fh(yi ~ as.factor(MajorArea), vardir = "var", combined_data = milk,
                                   domains = "dom", method = "amrl"))
  # emdi's optimize() tolerance (check S2)
  expect_equal(fit$sigma2_v, ref$model$variance, tolerance = 1e-3)
  expect_gt(fit$sigma2_v, scfh_sae(yi ~ as.factor(MajorArea), "var", data = milk, W = W, K = 1,
                                   method = "REML", mse = "none")$sigma2_v)
})

test_that("a recovered partition gives one Fay-Herriot fit per cluster", {
  skip_if_not_installed("sae")
  sim <- sim_regimes(b0 = c(0, 12), b1 = c(1, -1), seed = 2)
  d <- sim$data
  fit <- scfh_sae(y ~ x, "vardir", data = d, W = sim$W, K = 2, phi = 0.5, method = "REML", mse = "none", seed = 1)
  expect_equal(.adjusted_rand(fit$cluster$cluster, d$regime), 1)
  expect_true(fit$convergence)
  expect_true(fit$admissible)
  for (k in 1:2) {
    idx <- fit$cluster$cluster == k
    ref <- sae::eblupFH(y ~ x, vardir, data = d[idx, ], PRECISION = 1e-10)
    expect_equal(fit$sigma2_v[k], ref$fit$refvar, tolerance = 1e-6)
    expect_equal(fit$df_est$est[idx], unname(drop(ref$eblup)), tolerance = 1e-6)
    expect_equal(fit$coefficients[paste0("[", k, "] ", c("(Intercept)", "x")), "Estimate"],
                 unname(ref$fit$estcoef$beta), tolerance = 1e-6)
  }
  # the objective is the log-likelihood plus phi times the within-cluster edges
  cl <- fit$cluster$cluster
  edges <- sum(sim$W * outer(cl, cl, "==")) / 2
  expect_equal(fit$objective, fit$loglik + 0.5 * edges)
  expect_equal(sum(fit$size), nrow(d))
  expect_identical(fit$size[1] >= fit$size[2], TRUE)
})

test_that("the spatial penalty recovers regimes and smooths the partition", {
  sim <- sim_regimes(seed = 3)
  d <- sim$data
  cut_edges <- function(fit) {
    cl <- fit$cluster$cluster
    sum(sim$W * outer(cl, cl, "!=")) / 2
  }
  # moderate penalties: a very large one can lock the greedy label step into
  # the initial (covariate-based) partition
  fits <- lapply(c(0, 0.5, 2), function(phi) {
    scfh_sae(y ~ x, "vardir", data = d, W = sim$W, K = 2, phi = phi, mse = "none", seed = 1)
  })
  edges <- vapply(fits, cut_edges, numeric(1))
  expect_true(all(diff(edges) <= 0))
  expect_gte(.adjusted_rand(fits[[2]]$cluster$cluster, d$regime), 0.6)
  # well separated regimes are recovered
  clear <- sim_regimes(b0 = c(0, 10), b1 = c(2, -2), seed = 3)
  fit <- scfh_sae(y ~ x, "vardir", data = clear$data, W = clear$W, K = 2, phi = 0.5, mse = "none", seed = 1)
  expect_gte(.adjusted_rand(fit$cluster$cluster, clear$data$regime), 0.95)
})

test_that("results are reproducible with a seed and the session stream is restored", {
  sim <- sim_regimes(seed = 4)
  withr::local_seed(10)
  before <- runif(1)
  withr::local_seed(10)
  a <- scfh_sae(y ~ x, "vardir", data = sim$data, W = sim$W, K = 2, B = 5, seed = 7)
  expect_identical(runif(1), before)
  b <- scfh_sae(y ~ x, "vardir", data = sim$data, W = sim$W, K = 2, B = 5, seed = 7)
  expect_identical(a$df_est, b$df_est)
  expect_identical(a$boot_cluster, b$boot_cluster)
})

test_that("unsampled areas take their neighbours' cluster and the synthetic estimate", {
  sim <- sim_regimes(b0 = c(0, 10), b1 = c(2, -2), seed = 5)
  d <- sim$data
  d$y[c(12, 89)] <- NA
  W <- sim$W
  island <- 45
  W[island, ] <- W[, island] <- 0
  d$y[island] <- NA
  fit <- scfh_sae(y ~ x, "vardir", data = d, W = W, K = 2, phi = 0.5, B = 10, seed = 1)
  cl <- fit$cluster$cluster
  for (i in c(12, 89)) {
    expect_identical(cl[i], which.max(tabulate(cl[W[i, ] > 0], nbins = 2)))
  }
  expect_false(anyNA(cl))
  de <- fit$df_est
  expect_equal(de$est[!de$sampled], de$synthetic[!de$sampled])
  expect_true(all(is.na(de$gamma[!de$sampled])))
  expect_false(anyNA(de$mse))
  # row-standardized weights give the same neighbours
  Ws <- W / pmax(rowSums(W), 1)
  fit2 <- scfh_sae(y ~ x, "vardir", data = d, W = Ws, K = 2, phi = 0.5, mse = "none", seed = 1)
  expect_identical(fit2$cluster$cluster, cl)
})

test_that("bootstrap labels are aligned and the replicates feed benchmarking", {
  sim <- sim_regimes(b0 = c(0, 10), b1 = c(2, -2), seed = 6)
  fit <- scfh_sae(y ~ x, "vardir", data = sim$data, W = sim$W, K = 2, B = 10, seed = 1)
  expect_identical(dim(fit$boot_cluster), c(100L, 10L))
  expect_gt(mean(fit$boot_cluster == fit$cluster$cluster), 0.9)
  expect_equal(fit$df_est$mse, rowMeans(fit$boot_errors^2))
  expect_identical(fit$B_success, 10L)
  bm <- benchmark_flexsae(fit, target = 5)
  expect_equal(mean(bm$benchmarked), 5)
  expect_false(anyNA(bm$mse_benchmarked))
})

test_that(".align_labels() and .adjusted_rand() follow their definitions", {
  ref <- c(1, 1, 1, 2, 2, 2)
  expect_identical(.align_labels(c(2L, 2L, 2L, 1L, 1L, 1L), ref, 2), c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_identical(.align_labels(c(2L, 2L, 1L, 1L, 1L, 1L), ref, 2), c(1L, 1L, 2L, 2L, 2L, 2L))
  # Hubert and Arabie (1985): (2 - 1.2) / (4.5 - 1.2)
  expect_equal(.adjusted_rand(c(1, 1, 1, 2, 2, 2), c(1, 1, 2, 2, 3, 3)), 0.8 / 3.3)
  expect_equal(.adjusted_rand(c(1, 1, 2, 2), c(2, 2, 1, 1)), 1)
})

test_that("grid selection follows the two-step rule", {
  sim <- sim_regimes(b0 = c(0, 10), b1 = c(2, -2), seed = 7)
  fit <- scfh_sae(y ~ x, "vardir", data = sim$data, W = sim$W, K = 2, mse = "none", seed = 1)
  sel <- select_flexsae(fit, K = 1:3, phi = c(0, 0.25, 0.5, 0.75, 1), trace = FALSE)
  g <- sel$grid
  expect_identical(nrow(g), 11L)
  expect_s3_class(sel$model, "flexsae_scfh")
  expect_identical(sel$model$K, as.integer(sel$selected[["K"]]))
  # step 1 by hand
  band <- c(0.25, 0.5, 0.75, 1)
  pick <- vapply(band, function(p) {
    r <- g[(g$phi == p | g$K == 1) & g$admissible, ]
    r$K[which.min(r$criterion)]
  }, numeric(1))
  freq <- table(pick)
  expect_equal(sel$selected[["K"]], min(as.numeric(names(freq)[freq == max(freq)])))
  expect_equal(g$criterion, -2 * g$loglik + 3 * g$K * log(100) + ifelse(g$admissible, 0, NA))
  expect_output(print(sel), "Selected: K =")
  expect_error(select_flexsae(fit, method = "step"), "grid")
  expect_error(select_flexsae(fit, band = c(2, 3)), "band")
})

test_that("print, summary and plots show the clusters", {
  sim <- sim_regimes(seed = 8)
  fit <- scfh_sae(y ~ x, "vardir", data = sim$data, W = sim$W, K = 2, B = 5, seed = 1)
  expect_output(print(fit), "Clusters: K = 2")
  expect_output(print(summary(fit)), "sigma2_v")
  expect_equal(residuals(fit, "pearson"),
               stats::setNames((fit$df_est$direct - fit$df_est$synthetic) /
                                 sqrt(fit$sigma2_v[fit$cluster$cluster] + fit$df_est$vardir), fit$df_est$domain))
  for (type in c("clusters", "estimates", "residuals", "caterpillar", "rse")) {
    expect_silent(ggplot2::ggplot_build(autoplot(fit, type = type)))
  }
  expect_error(partial_flexsae(fit), "clusters")
  expect_error(autoplot(pspline_sae(y ~ x, "vardir", data = sim$data, mse = "none"), type = "clusters"),
               "scfh_sae")
})

test_that("invalid input gives clear errors", {
  sim <- sim_regimes(seed = 9)
  d <- sim$data
  W <- sim$W
  expect_error(scfh_sae(y ~ x, "vardir", data = d, K = 2), "required")
  expect_error(scfh_sae(y ~ x, "vardir", data = d, W = W[-1, ], K = 2), "100 x 100")
  W2 <- W
  W2[1, 2] <- -1
  expect_error(scfh_sae(y ~ x, "vardir", data = d, W = W2, K = 2), "negative")
  W2[1, 2] <- NA
  expect_error(scfh_sae(y ~ x, "vardir", data = d, W = W2, K = 2), "missing")
  expect_error(scfh_sae(y ~ x, "vardir", data = d, W = W, K = 1.5), "whole number")
  expect_error(scfh_sae(y ~ x, "vardir", data = d, W = W, phi = -1), "non-negative")
  expect_error(scfh_sae(y ~ s(x), "vardir", data = d, W = W), "smooth")
  expect_error(scfh_sae(y ~ x, "vardir", data = d, W = W, K = 40), "sampled areas")
})

test_that("a factor that is constant within a cluster is dropped from that cluster's fit", {
  sim <- sim_regimes(b0 = c(0, 6), b1 = c(1, -1), seed = 10)
  d <- sim$data
  d$f <- factor(rep(c("a", "b"), 50))
  # the initial k-means splits on the dummy, so each cluster sees one level
  expect_no_warning(fit <- scfh_sae(y ~ x + f, "vardir", data = d, W = sim$W, K = 2, mse = "none", seed = 1))
  expect_true(fit$admissible)
  expect_false(anyNA(fit$df_est$est))
  cf <- fit$coefficients
  expect_true(all(c("[1] fb", "[2] fb") %in% rownames(cf)))
  # an aliased column has an NA coefficient and contributes nothing
  for (k in 1:2) {
    idx <- fit$cluster$cluster == k
    X <- model.matrix(~ x + f, d)[idx, ]
    b <- cf[paste0("[", k, "] ", colnames(X)), "Estimate"]
    expect_equal(fit$df_est$synthetic[idx], unname(drop(X %*% ifelse(is.na(b), 0, b))))
    if (length(unique(d$f[idx])) == 1) expect_true(is.na(b[3]))
  }
  # the aliased path itself: a column constant in the cluster
  X <- cbind("(Intercept)" = 1, x = d$x, fb = 1)
  a <- .fit_fh_aliased(d$y, X, d$vardir, "REML")
  ref <- .fit_fh_lm(d$y, X[, 1:2], d$vardir, "REML")
  expect_identical(unname(a$beta), c(unname(ref$beta), 0))
  expect_true(is.na(a$coefficients["fb", "Estimate"]))
  expect_equal(a$sigma2_v, ref$sigma2_v)
})
