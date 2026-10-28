# Fay-Herriot fit through mgcv, shared by the area-level models.
#
# y_d = m(x_d) + v_d + e_d,  v_d ~ N(0, sigma2_v),  e_d ~ N(0, vardir_d) known.
# Known sampling variances: weights 1 / vardir with the scale fixed at 1; the
# area effect is an `s(area, bs = "re")` smooth, whose variance mgcv estimates
# as scale / sp by REML or ML (Wood, 2011). With a linear m() this reproduces
# sae::eblupFH (technical-report/DISCREPANCIES.md, check S1).

#' Fit the Fay-Herriot model with mgcv on the sampled areas
#'
#' @param formula Model formula from `.fh_formula()` (internal response and
#'   area term already in place).
#' @param fit_data Data of the sampled areas, with the internal response,
#'   area factor and weight columns.
#' @param method "REML" or "ML".
#' @return List with the `gam` fit, `sigma2_v` and `converged`.
#' @noRd
.fit_fh_gam <- function(formula, fit_data, method) {
  fit <- tryCatch(
    mgcv::gam(formula, data = fit_data, weights = .flexsae_w, scale = 1, method = method),
    error = function(e) {
      cli::cli_abort(c(
        "The model could not be fitted.",
        "x" = conditionMessage(e),
        "i" = "With few areas, lower the basis dimension of the smooths, e.g. {.code s(x, k = 5)}."
      ), call = NULL)
    }
  )
  sp_name <- paste0("s(", .area_var, ")")
  list(
    fit = fit,
    sigma2_v = unname(fit$scale / fit$sp[[sp_name]]),
    converged = identical(fit$outer.info$conv, "full convergence")
  )
}

#' Synthetic part and area effects of a Fay-Herriot gam fit for all areas
#'
#' @param fh Result of `.fit_fh_gam()`.
#' @param newdata Data of all areas (sampled and unsampled).
#' @param sampled Logical vector over the rows of `newdata`.
#' @return List with `synthetic` (m-hat for every area) and `random_effect`
#'   (v-hat, 0 for unsampled areas).
#' @noRd
.predict_fh_gam <- function(fh, newdata, sampled) {
  fit <- fh$fit
  lev <- levels(fit$model[[.area_var]])
  # the area term is excluded below; unsampled rows only need a valid level
  area <- as.character(seq_len(nrow(newdata)))
  area[!sampled] <- lev[1]
  newdata[[.area_var]] <- factor(area, levels = lev)
  synthetic <- tryCatch(
    as.numeric(mgcv::predict.gam(fit, newdata = newdata, exclude = paste0("s(", .area_var, ")"))),
    error = function(e) {
      cli::cli_abort(c(
        "Prediction for the areas failed.",
        "x" = conditionMessage(e),
        "i" = "Unsampled areas need covariate values (and factor levels) seen in the sampled areas."
      ), call = NULL)
    }
  )
  coefs <- stats::coef(fit)
  v <- coefs[startsWith(names(coefs), paste0("s(", .area_var, ")"))]
  random_effect <- numeric(nrow(newdata))
  random_effect[as.integer(lev)] <- unname(v)
  list(synthetic = synthetic, random_effect = random_effect)
}
