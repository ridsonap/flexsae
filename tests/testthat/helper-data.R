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
