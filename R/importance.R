#' Variable Importance of a Machine Learning Fay-Herriot Model
#'
#' Importance of the covariates in the final forest of a [merf_sae()] fit or
#' the final booster of a [boost_sae()] fit.
#'
#' @param object A model from [merf_sae()] or [boost_sae()].
#' @param type For [merf_sae()]: `"permutation"` (default), the increase of
#'   the out-of-bag mean squared error when the values of a covariate are
#'   permuted (Breiman, 2001), or `"impurity"`, the total decrease of node
#'   impurity. For [boost_sae()]: `"gain"` (the default and only choice), the
#'   total gain of the splits on a covariate, summed over the one-hot columns
#'   of a factor and normalized to sum to one (`xgboost::xgb.importance()`).
#'
#' @details The forest is refitted to the stored response of the final forest
#' (`rf_target`) with the same seed and settings, so its trees and predictions
#' equal those of the model; ranger then computes the importance. Permutation
#' importance can favour correlated covariates and impurity importance
#' covariates with many split points (Strobl et al., 2008).
#'
#' @return A data frame with columns `variable` and `importance`, sorted by
#'   decreasing importance, with attribute `type`.
#'
#' @references
#' Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5-32.
#' \doi{10.1023/A:1010933404324}
#'
#' Strobl, C., Boulesteix, A.-L., Kneib, T., Augustin, T. and Zeileis, A.
#' (2008). Conditional variable importance for random forests. *BMC
#' Bioinformatics*, 9, 307. \doi{10.1186/1471-2105-9-307}
#'
#' @export
#' @examples
#' set.seed(1)
#' D <- 80
#' dat <- data.frame(x1 = runif(D, -2, 2), x2 = rnorm(D), x3 = rnorm(D))
#' dat$vardir <- runif(D, 0.2, 0.6)
#' dat$y <- 2 + sin(1.5 * dat$x1) + 0.5 * dat$x2 + rnorm(D, sd = 0.5) +
#'   rnorm(D, sd = sqrt(dat$vardir))
#' fit <- merf_sae(y ~ x1 + x2 + x3, vardir = "vardir", data = dat,
#'                 num.trees = 200, mse = "none", seed = 1)
#' importance_flexsae(fit)
importance_flexsae <- function(object, type = c("permutation", "impurity", "gain")) {
  if (inherits(object, "flexsae_boost")) {
    if (!missing(type) && !identical(type, "gain")) {
      cli::cli_abort("{.fn boost_sae} models have {.code type = \"gain\"} importance only.", call = NULL)
    }
    return(.boost_importance(object))
  }
  type <- if (missing(type)) "permutation" else match.arg(type)
  if (type == "gain" && inherits(object, "flexsae_merf")) {
    cli::cli_abort("{.code type = \"gain\"} is for {.fn boost_sae} models.", call = NULL)
  }
  if (!inherits(object, "flexsae_merf")) {
    cli::cli_abort(c(
      "{.fn importance_flexsae} needs a {.fn merf_sae} or {.fn boost_sae} model.",
      "i" = "For {.fn pspline_sae} models see {.code summary(fit$fit)} and {.fn partial_flexsae}."
    ), call = NULL)
  }
  rf <- object$rf
  x <- .rf_covariates(object)[object$df_est$sampled, , drop = FALSE]
  fit <- ranger::ranger(
    x = x, y = object$rf_target, num.trees = rf$num.trees, mtry = rf$mtry,
    min.node.size = rf$min.node.size, seed = object$rf_seed, num.threads = 1,
    importance = type
  )
  imp <- fit$variable.importance
  out <- data.frame(variable = names(imp), importance = unname(imp), stringsAsFactors = FALSE)
  out <- out[order(-out$importance), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "type") <- type
  out
}

#' Covariates of the forest of a merf_sae() fit, all areas
#' @noRd
.rf_covariates <- function(object) {
  x <- object$data[, all.vars(object$formula[[3]]), drop = FALSE]
  x[] <- lapply(x, function(v) if (is.character(v)) factor(v) else v)
  x
}

#' Gain importance of a boost_sae() booster, by covariate
#' @noRd
.boost_importance <- function(object) {
  ba <- object$boost_args
  gain <- stats::setNames(numeric(length(ba$columns)), ba$columns)
  imp <- xgboost::xgb.importance(model = object$booster)
  if (nrow(imp) > 0) gain[imp$Feature] <- imp$Gain
  by_var <- tapply(gain, factor(ba$variables[ba$assign], levels = ba$variables), sum)
  out <- data.frame(variable = names(by_var), importance = as.numeric(by_var), stringsAsFactors = FALSE)
  out <- out[order(-out$importance), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "type") <- "gain"
  out
}

