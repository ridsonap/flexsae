# Data-generating processes of the simulation study (PLAN.md section 3.2).
# Area-level model: theta_d = f(x_d) + v_d, y_d = theta_d + e_d, e_d ~ N(0, psi_d).
# Covariates are drawn once per scenario (fixed design); v and e are redrawn in
# every replicate. All mean functions are our own designs except A5, the
# non-linear mean function of Bukhari et al. (2026, eq. 20), with covariate
# distributions chosen by us (the paper does not give them in the text read).

mean_functions <- list(
  B = function(x, lambda = 1) 1 + x$x1 + 0.5 * x$x2,
  A1 = function(x, lambda = 1) 1 + 0.8 * x$x1^2 + 0.5 * x$x2,
  A2 = function(x, lambda = 1) 1 + lambda * 1.5 * sin(1.5 * x$x1) + (1 - lambda) * x$x1 + 0.5 * x$x2,
  A3 = function(x, lambda = 1) 1 + 1.5 * (x$x1 > 0) + 0.5 * x$x2,
  A4 = function(x, lambda = 1) 1 + 0.8 * x$x1 * x$x2 + 0.5 * x$x2,
  A5 = function(x, lambda = 1) x$x1^2 + 2 * x$x2^3 + 3 * log(abs(x$x3) + 1) + 4 * log(abs(x$x4)),
  C = function(x, lambda = 1) 1 + 1.5 * sin(1.5 * x$x1) + 0.8 * x$x2^2 + 0.5 * x$x3
)

# covariates of the fitted models, by mean function
model_covariates <- function(shape) {
  switch(shape,
    A5 = paste0("x", 1:4),
    C = paste0("x", 1:8),
    c("x1", "x2")
  )
}

scenario <- function(id, shape, D = 100, ratio = 1, psi = "uniform", v_dist = "normal",
                     unsampled = 0, lambda = 1, sigma2_v = 1) {
  list(id = id, shape = shape, D = D, ratio = ratio, psi = psi, v_dist = v_dist,
       unsampled = unsampled, lambda = lambda, sigma2_v = sigma2_v)
}

point_scenarios <- function() {
  base <- list(
    scenario("B", "B"), scenario("A1", "A1"), scenario("A2", "A2"), scenario("A3", "A3"),
    scenario("A4", "A4"), scenario("A5", "A5"), scenario("C", "C"),
    scenario("A2_lambda0.5", "A2", lambda = 0.5)
  )
  vary <- function(s) {
    list(
      scenario(paste0(s, "_D30"), s, D = 30), scenario(paste0(s, "_D50"), s, D = 50),
      scenario(paste0(s, "_D300"), s, D = 300),
      scenario(paste0(s, "_ratio0.25"), s, ratio = 0.25), scenario(paste0(s, "_ratio4"), s, ratio = 4),
      scenario(paste0(s, "_psi5"), s, psi = "groups"), scenario(paste0(s, "_t3"), s, v_dist = "t3"),
      scenario(paste0(s, "_unsampled"), s, unsampled = 0.1)
    )
  }
  c(base, vary("B"), vary("A2"))
}

mse_scenarios <- function() {
  list(scenario("B", "B"), scenario("A1", "A1"), scenario("A2", "A2"), scenario("C", "C"),
       scenario("A2_D50", "A2", D = 50))
}

# fixed part of a scenario: covariates, f(x), sampling variances
make_design <- function(sc, seed) {
  withr::with_seed(seed, {
    D <- sc$D
    if (sc$shape == "A5") {
      x <- as.data.frame(matrix(runif(D * 4, 0.5, 2), D, 4, dimnames = list(NULL, paste0("x", 1:4))))
    } else if (sc$shape == "C") {
      x <- as.data.frame(matrix(runif(D * 8, -2, 2), D, 8, dimnames = list(NULL, paste0("x", 1:8))))
      x$x4 <- 0.7 * x$x1 + sqrt(1 - 0.7^2) * x$x4 # noise correlated with x1 (r about 0.7)
      x$x5 <- 0.7 * x$x3 + sqrt(1 - 0.7^2) * x$x5 # noise correlated with x3
    } else {
      x <- data.frame(x1 = runif(D, -2, 2), x2 = runif(D, -2, 2))
    }
    f <- mean_functions[[sc$shape]](x, sc$lambda)
    # sampling variances with mean psi_bar = sigma2_v / ratio
    psi_bar <- sc$sigma2_v / sc$ratio
    psi <- if (sc$psi == "groups") {
      psi_bar * rep(c(0.1, 0.3, 0.8, 1.3, 2.5), length.out = D)
    } else {
      psi_bar * runif(D, 0.5, 1.5)
    }
    sampled <- rep(TRUE, D)
    if (sc$unsampled > 0) sampled[seq_len(round(sc$unsampled * D))] <- FALSE
    list(x = x, f = f, psi = psi, sampled = sampled)
  })
}

# one replicate: theta and direct estimates
draw_replicate <- function(sc, design) {
  D <- sc$D
  v <- if (sc$v_dist == "t3") sqrt(sc$sigma2_v / 3) * rt(D, 3) else rnorm(D, 0, sqrt(sc$sigma2_v))
  theta <- design$f + v
  y <- theta + rnorm(D, 0, sqrt(design$psi))
  y[!design$sampled] <- NA
  dat <- design$x
  dat$y <- y
  dat$psi <- design$psi
  list(data = dat, theta = theta)
}
