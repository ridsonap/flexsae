#' @importFrom ggplot2 autoplot .data
#' @export
ggplot2::autoplot

#' Plot a flexsae Model
#'
#' `ggplot2` plots of the estimates, their precision and the model diagnostics.
#'
#' @param object A `"flexsae"` model.
#' @param type Plot type:
#'   * `"estimates"` (default): direct and model estimates by area, ordered by
#'     the sampling variance, with the confidence intervals of the model;
#'   * `"scatter"`: model against direct estimates, with the identity line;
#'   * `"rse"`: relative standard errors (\%, log scale) of the direct and model
#'     estimates against the sampling variance, with reference lines at 25 and
#'     50\%;
#'   * `"residuals"`: standardized residuals
#'     \eqn{(y_d - \hat m_d)/\sqrt{\hat\sigma^2_v + \psi_d}} against \eqn{\hat m_d};
#'   * `"qq"`: normal quantile plot of the standardized residuals;
#'   * `"caterpillar"`: \eqn{\hat v_d} with intervals
#'     \eqn{\pm z\sqrt{\gamma_d\psi_d}} (the conditional standard deviation of
#'     \eqn{v_d} with known parameters);
#'   * `"partial"`: partial dependence of the fixed part ([partial_flexsae()]);
#'   * `"importance"`: variable importance ([importance_flexsae()], `merf_sae()`
#'     models);
#'   * `"clusters"`: coefficients by cluster with intervals
#'     \eqn{\pm z\,\mathrm{SE}} ([scfh_sae()] models).
#' @param ... Passed to [partial_flexsae()] (`type = "partial"`) or
#'   [importance_flexsae()] (`type = "importance"`).
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' set.seed(1)
#' D <- 60
#' dat <- data.frame(x1 = runif(D, -2, 2), x2 = rnorm(D))
#' dat$vardir <- runif(D, 0.2, 0.6)
#' dat$y <- 2 + sin(1.5 * dat$x1) + 0.5 * dat$x2 + rnorm(D, sd = 0.5) +
#'   rnorm(D, sd = sqrt(dat$vardir))
#' fit <- pspline_sae(y ~ s(x1, k = 6) + x2, vardir = "vardir", data = dat,
#'                    B = 20, seed = 1)
#' autoplot(fit)
#' autoplot(fit, type = "partial")
autoplot.flexsae <- function(object, type = c("estimates", "scatter", "rse", "residuals", "qq",
                                              "caterpillar", "partial", "importance", "clusters"), ...) {
  type <- match.arg(type)
  df <- object$df_est
  z <- stats::qnorm(1 - (1 - object$conf_level) / 2)
  theme <- ggplot2::theme_minimal(base_size = 11)
  col_direct <- "#9A9A9A"
  col_model <- "#1D6A5C"
  s <- df$sampled

  if (type == "clusters" && !inherits(object, "flexsae_scfh")) {
    cli::cli_abort("{.code type = \"clusters\"} needs a {.fn scfh_sae} model.", call = NULL)
  }
  if (type == "rse" && all(is.na(df$mse))) {
    cli::cli_abort("The model has no MSE; fit it with {.code mse = \"bootstrap\"}.", call = NULL)
  }

  p <- switch(type,
    estimates = {
      ord <- order(df$vardir, na.last = TRUE)
      long <- rbind(
        data.frame(rank = seq_along(ord), value = df$direct[ord], lower = NA_real_, upper = NA_real_,
                   estimate = "Direct"),
        data.frame(rank = seq_along(ord), value = df$est[ord], lower = df$ci_lower[ord],
                   upper = df$ci_upper[ord], estimate = "Model")
      )
      g <- ggplot2::ggplot(long, ggplot2::aes(x = .data$rank, y = .data$value, colour = .data$estimate))
      if (!all(is.na(long$lower))) {
        g <- g + ggplot2::geom_errorbar(
          data = long[long$estimate == "Model", ],
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper), width = 0, alpha = 0.5
        )
      }
      g + ggplot2::geom_point(size = 1.6, na.rm = TRUE) +
        ggplot2::scale_colour_manual(values = c(Direct = col_direct, Model = col_model)) +
        ggplot2::labs(x = "Area (by increasing sampling variance)", y = "Estimate", colour = NULL)
    },
    scatter = {
      d <- df[s, ]
      ggplot2::ggplot(d, ggplot2::aes(x = .data$direct, y = .data$est)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50") +
        ggplot2::geom_point(colour = col_model, alpha = 0.8) +
        ggplot2::labs(x = "Direct estimate", y = "Model estimate")
    },
    rse = {
      rse_direct <- 100 * sqrt(df$vardir) / abs(df$direct)
      long <- rbind(
        data.frame(vardir = df$vardir, rse = rse_direct, estimate = "Direct"),
        data.frame(vardir = df$vardir, rse = df$rse, estimate = "Model")
      )
      long <- long[is.finite(long$rse) & is.finite(long$vardir), ]
      ggplot2::ggplot(long, ggplot2::aes(x = .data$vardir, y = .data$rse, colour = .data$estimate)) +
        ggplot2::geom_hline(yintercept = c(25, 50), linetype = "dashed", colour = "grey60") +
        ggplot2::geom_point(alpha = 0.8) +
        ggplot2::scale_colour_manual(values = c(Direct = col_direct, Model = col_model)) +
        ggplot2::scale_y_log10() +
        ggplot2::labs(x = "Sampling variance", y = "Relative standard error (%, log scale)", colour = NULL)
    },
    residuals = {
      d <- data.frame(fitted = df$synthetic[s], resid = stats::residuals(object, type = "pearson")[s])
      ggplot2::ggplot(d, ggplot2::aes(x = .data$fitted, y = .data$resid)) +
        ggplot2::geom_hline(yintercept = c(-2, 0, 2), linetype = c("dashed", "solid", "dashed"),
                            colour = "grey60") +
        ggplot2::geom_point(colour = col_model, alpha = 0.8) +
        ggplot2::labs(x = "Fixed part", y = "Standardized residual")
    },
    qq = {
      d <- data.frame(resid = stats::residuals(object, type = "pearson")[s])
      ggplot2::ggplot(d, ggplot2::aes(sample = .data$resid)) +
        ggplot2::stat_qq_line(colour = "grey50") +
        ggplot2::stat_qq(colour = col_model) +
        ggplot2::labs(x = "Normal quantile", y = "Standardized residual")
    },
    caterpillar = {
      d <- data.frame(v = df$random_effect[s], sd = sqrt(df$gamma[s] * df$vardir[s]))
      d <- d[order(d$v), ]
      d$rank <- seq_len(nrow(d))
      ggplot2::ggplot(d, ggplot2::aes(x = .data$rank, y = .data$v)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$v - z * .data$sd, ymax = .data$v + z * .data$sd),
                               width = 0, colour = col_model, alpha = 0.5) +
        ggplot2::geom_point(colour = col_model, size = 1.4) +
        ggplot2::labs(x = "Area (by estimated effect)", y = "Area random effect")
    },
    partial = {
      pd <- partial_flexsae(object, ...)
      num <- pd[!is.na(pd$x), ]
      if (nrow(num) == 0) {
        cli::cli_abort("{.code type = \"partial\"} plots numeric covariates; none were given.", call = NULL)
      }
      g <- ggplot2::ggplot(num, ggplot2::aes(x = .data$x, y = .data$fit))
      ice <- attr(pd, "ice")
      if (!is.null(ice)) {
        ice <- ice[!is.na(ice$x), ]
        g <- g + ggplot2::geom_line(data = ice, ggplot2::aes(group = .data$domain), colour = "grey80",
                                    linewidth = 0.3)
      }
      if (!all(is.na(num$se))) {
        g <- g + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                                      fill = col_model, alpha = 0.2)
      }
      g + ggplot2::geom_line(colour = col_model, linewidth = 0.9) +
        ggplot2::facet_wrap(~variable, scales = "free_x") +
        ggplot2::labs(x = NULL, y = "Partial dependence of the fixed part")
    },
    importance = {
      imp <- importance_flexsae(object, ...)
      imp$variable <- factor(imp$variable, levels = rev(imp$variable))
      ggplot2::ggplot(imp, ggplot2::aes(x = .data$importance, y = .data$variable)) +
        ggplot2::geom_col(fill = col_model) +
        ggplot2::labs(x = paste0("Importance (", attr(imp, "type"), ")"), y = NULL)
    },
    clusters = {
      cf <- object$coefficients
      d <- data.frame(
        cluster = factor(sub("^\\[(\\d+)\\] .*", "\\1", rownames(cf)), levels = seq_len(object$K)),
        term = sub("^\\[\\d+\\] ", "", rownames(cf)),
        estimate = cf$Estimate,
        se = cf$`Std. Error`
      )
      ggplot2::ggplot(d, ggplot2::aes(x = .data$estimate, y = .data$cluster)) +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$estimate - z * .data$se,
                                            xmax = .data$estimate + z * .data$se),
                               width = 0, orientation = "y", colour = col_model, alpha = 0.6) +
        ggplot2::geom_point(colour = col_model, size = 1.8) +
        ggplot2::facet_wrap(~term, scales = "free_x") +
        ggplot2::labs(x = "Coefficient", y = "Cluster")
    }
  )
  p + theme
}

#' Plot Benchmarked Estimates
#'
#' @param object A `"flexsae_benchmark"` object from [benchmark_flexsae()].
#' @param type `"estimates"` (default): benchmarked against original
#'   estimates; `"mse"`: MSE after against before benchmarking.
#' @param ... Not used.
#' @return A `ggplot` object.
#' @export
autoplot.flexsae_benchmark <- function(object, type = c("estimates", "mse"), ...) {
  type <- match.arg(type)
  d <- as.data.frame(object)
  many_groups <- length(unique(d$group)) > 1
  if (type == "mse" && all(is.na(d$mse_benchmarked))) {
    cli::cli_abort("The benchmarked estimates have no MSE.", call = NULL)
  }
  aes_xy <- if (type == "estimates") {
    ggplot2::aes(x = .data$original, y = .data$benchmarked)
  } else {
    ggplot2::aes(x = .data$mse, y = .data$mse_benchmarked)
  }
  g <- ggplot2::ggplot(d, aes_xy) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50")
  g <- if (many_groups) {
    g + ggplot2::geom_point(ggplot2::aes(colour = .data$group), alpha = 0.85) +
      ggplot2::labs(colour = "Group")
  } else {
    g + ggplot2::geom_point(colour = "#1D6A5C", alpha = 0.85)
  }
  labs <- if (type == "estimates") {
    ggplot2::labs(x = "Original estimate", y = "Benchmarked estimate")
  } else {
    ggplot2::labs(x = "MSE before benchmarking", y = "MSE after benchmarking")
  }
  g <- g + labs + ggplot2::theme_minimal(base_size = 11)
  # a legend with many groups (e.g. 34 provinces) is not readable
  if (length(unique(d$group)) > 10) g <- g + ggplot2::theme(legend.position = "none")
  g
}
