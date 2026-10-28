# Area-level data with a non-linear covariate effect:
# theta_d = 2 + sin(1.5 x1) + 0.5 x2 + v_d,  y_d = theta_d + e_d.
sim_nonlinear_area <- function(D = 80, sigma_v = 0.5, seed = 1) {
  withr::with_seed(seed, {
    dat <- data.frame(area = paste0("A", seq_len(D)), x1 = runif(D, -2, 2), x2 = rnorm(D))
    dat$vardir <- runif(D, 0.2, 0.6)
    dat$theta <- 2 + sin(1.5 * dat$x1) + 0.5 * dat$x2 + rnorm(D, sd = sigma_v)
    dat$y <- dat$theta + rnorm(D, sd = sqrt(dat$vardir))
    dat
  })
}

milk_data <- function() {
  data("milk", package = "sae", envir = environment())
  milk$var <- milk$SD^2
  milk
}

# Two spatial regimes on an n_side x n_side lattice (left and right halves):
# theta_d = b0_k + b1_k x_d + u_d, rook adjacency in W.
sim_regimes <- function(n_side = 10, b0 = c(2, 6), b1 = c(2, -1), sigma_u = 0.5, seed = 1) {
  withr::with_seed(seed, {
    g <- expand.grid(col = seq_len(n_side), row = seq_len(n_side))
    D <- nrow(g)
    regime <- ifelse(g$col <= n_side / 2, 1L, 2L)
    dat <- data.frame(area = paste0("A", seq_len(D)), x = rnorm(D, sd = 2), regime = regime)
    dat$vardir <- runif(D, 0.3, 0.6)
    dat$theta <- b0[regime] + b1[regime] * dat$x + rnorm(D, sd = sigma_u)
    dat$y <- dat$theta + rnorm(D, sd = sqrt(dat$vardir))
    list(data = dat, W = 1 * (as.matrix(stats::dist(g)) == 1))
  })
}
