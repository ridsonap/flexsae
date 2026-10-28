#' Covariate Selection for flexsae Models
#'
#' Selects the covariates of a fitted model and refits it with them.
#'
#' @param object A model from [pspline_sae()] or [merf_sae()].
#' @param method `"step"` (stepwise search on an information criterion;
#'   `pspline_sae()` models) or `"boruta"` (all-relevant selection with the
#'   Boruta algorithm; `merf_sae()` models). The default depends on the model.
#' @param criterion For `"step"`: `"BIC"` (default), `"AIC"` or `"KIC"`.
#' @param direction For `"step"`: `"backward"` (default), `"both"` or
#'   `"forward"`.
#' @param scope For `"step"`: optional one-sided formula of the largest model,
#'   e.g. `~ s(x1) + s(x2) + x3`; defaults to the terms of `object`. The
#'   search starts from the model of `object`.
#' @param trace For `"step"`: print the criterion of each step.
#' @param max_runs,p_value For `"boruta"`: passed to `Boruta::Boruta()` as
#'   `maxRuns` and `pValue`.
#' @param seed For `"boruta"`: optional seed; the random number stream of the
#'   session is restored afterwards.
#'
#' @details
#' **Stepwise (`pspline_sae()`).** Every candidate model is fitted by maximum
#' likelihood; the log-likelihood \eqn{\ell} is the marginal Gaussian
#' log-likelihood of the sampled direct estimates, with the area effects and
#' the penalized spline coefficients integrated out (the ML score of
#' [mgcv::gam()]). With \eqn{k} the number of unpenalized coefficients plus the
#' number of variance parameters (\eqn{\sigma^2_v} and one per smooth term) and
#' \eqn{m} the number of sampled areas, \eqn{\mathrm{AIC} = -2\ell + 2k},
#' \eqn{\mathrm{BIC} = -2\ell + k\log m} and \eqn{\mathrm{KIC} = -2\ell + 3k}
#' (Marhuenda, Morales and Pardo, 2014). For linear terms only, \eqn{\ell}, the
#' criteria and the search are those of `fastsae::step_sae()`. For smooth terms
#' the criterion counts a smoothing parameter as one parameter, a simple
#' marginal criterion; it is a heuristic whose behaviour is studied by
#' simulation (not in the paper above). REML cannot compare models with
#' different fixed effects, so the search uses ML; the selected model is
#' refitted with the method, MSE and self-benchmarking settings of `object`.
#'
#' **Boruta (`merf_sae()`).** The Boruta algorithm (Kursa and Rudnicki, 2010)
#' compares the permutation importance of each covariate in a random forest
#' with that of permuted copies ("shadows"), here with `Boruta::getImpRfZ()`
#' (ranger). The forest is fitted to the target of the first step of the model,
#' \eqn{y_d - \hat v_d}, over the sampled areas. Covariates confirmed as
#' relevant are kept (tentative ones are dropped); the model is refitted with
#' them. Permutation importance can favour correlated covariates (Strobl et
#' al., 2008).
#'
#' @return An object of class `"flexsae_selection"`: a list with `model` (the
#'   refitted model), `formula`, `selected` (kept terms or variables),
#'   `method`, and for `"step"` the `path` (data frame of steps) and
#'   `criterion`, for `"boruta"` the `decision` table.
#'
#' @references
#' Kursa, M. B. and Rudnicki, W. R. (2010). Feature selection with the Boruta
#' package. *Journal of Statistical Software*, 36(11), 1-13.
#' \doi{10.18637/jss.v036.i11}
#'
#' Marhuenda, Y., Morales, D. and Pardo, M. C. (2014). Information criteria
#' for Fay-Herriot model selection. *Computational Statistics & Data
#' Analysis*, 70, 268-280. \doi{10.1016/j.csda.2013.09.016}
#'
#' Strobl, C., Boulesteix, A.-L., Kneib, T., Augustin, T. and Zeileis, A.
#' (2008). Conditional variable importance for random forests. *BMC
#' Bioinformatics*, 9, 307. \doi{10.1186/1471-2105-9-307}
#'
#' @export
#' @examples
#' set.seed(1)
#' D <- 60
#' dat <- data.frame(x1 = runif(D, -2, 2), x2 = rnorm(D), x3 = rnorm(D))
#' dat$vardir <- runif(D, 0.2, 0.6)
#' dat$y <- 2 + sin(1.5 * dat$x1) + 0.5 * dat$x2 + rnorm(D, sd = 0.5) +
#'   rnorm(D, sd = sqrt(dat$vardir))
#' fit <- pspline_sae(y ~ s(x1, k = 6) + x2 + x3, vardir = "vardir", data = dat,
#'                    mse = "none")
#' sel <- select_flexsae(fit)
#' sel
select_flexsae <- function(
  object,
  method = NULL,
  criterion = c("BIC", "AIC", "KIC"),
  direction = c("backward", "both", "forward"),
  scope = NULL,
  trace = TRUE,
  max_runs = 100,
  p_value = 0.01,
  seed = NULL
) {
  if (!inherits(object, "flexsae")) {
    cli::cli_abort("{.arg object} must be a model from {.fn pspline_sae} or {.fn merf_sae}.", call = NULL)
  }
  is_merf <- identical(object$model, "merf_sae")
  method <- method %||% if (is_merf) "boruta" else "step"
  method <- match.arg(method, c("step", "boruta"))
  if (method == "step" && is_merf) {
    cli::cli_abort(c(
      "{.code method = \"step\"} needs a likelihood and is available for {.fn pspline_sae} models.",
      "i" = "Use {.code method = \"boruta\"} for {.fn merf_sae} models."
    ), call = NULL)
  }
  if (method == "boruta" && !is_merf) {
    cli::cli_abort("{.code method = \"boruta\"} is available for {.fn merf_sae} models.", call = NULL)
  }
  criterion <- match.arg(criterion)
  direction <- match.arg(direction)
  if (method == "step") {
    .select_step(object, criterion, direction, scope, trace)
  } else {
    .select_boruta(object, max_runs, p_value, seed)
  }
}

#' Refit a model with a new formula and settings taken from the object
#' @noRd
.refit_flexsae <- function(object, formula, method = object$method, mse = object$mse_method) {
  args <- list(
    formula = formula, vardir = object$df_est$vardir, data = object$data,
    method = method, mse = mse, seed = object$seed, conf_level = object$conf_level,
    bench_weight = object$bench$weight, bench_group = object$bench$group
  )
  if (mse == "bootstrap") args$B <- object$B
  if (identical(object$model, "merf_sae")) {
    ra <- object$rf_args
    if (!is.null(ra$mtry) && ra$mtry > length(all.vars(formula[[3]]))) ra$mtry <- NULL
    args <- c(args, list(
      variant = object$variant, algorithm = object$algorithm,
      num.trees = ra$num.trees, mtry = ra$mtry, min.node.size = ra$min.node.size
    ))
    fit <- do.call(merf_sae, args)
  } else {
    fit <- do.call(pspline_sae, args)
  }
  fit$df_est$domain <- object$df_est$domain
  fit
}

#' Marginal log-likelihood and number of parameters of an ML pspline fit
#' @noRd
.pspline_ic <- function(fit, criterion) {
  g <- fit$fit
  ll <- -as.numeric(g$gcv.ubre) # drop the attributes mgcv attaches
  n_unpen <- g$nsdf
  for (sm in g$smooth) {
    S <- Reduce(`+`, lapply(sm$S, function(s) s / max(abs(s))))
    n_unpen <- n_unpen + (sm$last.para - sm$first.para + 1) - qr(S)$rank
  }
  k <- n_unpen + length(g$sp)
  m <- sum(fit$df_est$sampled)
  switch(criterion,
    AIC = -2 * ll + 2 * k,
    BIC = -2 * ll + k * log(m),
    KIC = -2 * ll + 3 * k
  )
}

#' Stepwise search for pspline_sae()
#' @noRd
.select_step <- function(object, criterion, direction, scope, trace) {
  f0 <- object$formula
  resp <- f0[[2]]
  current <- attr(stats::terms(f0), "term.labels")
  upper <- if (is.null(scope)) current else attr(stats::terms(scope), "term.labels")
  missing_terms <- setdiff(current, upper)
  if (length(missing_terms) > 0) {
    cli::cli_abort("{.arg scope} must contain the terms of the model: {.val {missing_terms}} {?is/are} missing.", call = NULL)
  }
  intercept <- attr(stats::terms(f0), "intercept") == 1
  build <- function(terms) {
    rhs <- if (length(terms) == 0) "1" else paste(terms, collapse = " + ")
    if (!intercept) rhs <- paste(rhs, "- 1")
    f <- stats::as.formula(paste(deparse(resp), "~", rhs), env = environment(f0))
    f
  }
  score <- function(terms) {
    fit <- suppressWarnings(.refit_flexsae(object, build(terms), method = "ML", mse = "none"))
    .pspline_ic(fit, criterion)
  }

  cur_score <- score(current)
  if (trace) cli::cli_text("Start: {criterion} = {round(cur_score, 4)} ({deparse(build(current))})")
  path <- data.frame(step = "start", change = "", value = cur_score, stringsAsFactors = FALSE)
  repeat {
    cands <- list()
    if (direction %in% c("backward", "both")) {
      for (t in current) cands[[paste("-", t)]] <- setdiff(current, t)
    }
    if (direction %in% c("forward", "both")) {
      for (t in setdiff(upper, current)) cands[[paste("+", t)]] <- c(current, t)
    }
    if (length(cands) == 0) break
    vals <- vapply(cands, score, numeric(1))
    best <- which.min(vals)
    if (vals[best] >= cur_score - 1e-10) break
    current <- cands[[best]]
    cur_score <- vals[[best]]
    path <- rbind(path, data.frame(step = as.character(nrow(path)), change = names(cands)[best],
                                   value = cur_score, stringsAsFactors = FALSE))
    if (trace) cli::cli_text("Step: {names(cands)[best]}: {criterion} = {round(cur_score, 4)}")
  }
  names(path)[3] <- criterion
  final_formula <- build(current)
  model <- .refit_flexsae(object, final_formula)
  structure(
    list(model = model, formula = final_formula, selected = current, method = "step",
         path = path, criterion = criterion, direction = direction),
    class = "flexsae_selection"
  )
}

#' Boruta selection for merf_sae()
#' @noRd
.select_boruta <- function(object, max_runs, p_value, seed) {
  if (!requireNamespace("Boruta", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg Boruta} is needed for {.code method = \"boruta\"}.", call = NULL)
  }
  f0 <- object$formula
  vars <- all.vars(f0[[3]])
  if (length(vars) < 2) {
    cli::cli_abort("Selection needs at least two covariates.", call = NULL)
  }
  s <- object$df_est$sampled
  x <- object$data[s, vars, drop = FALSE]
  x[] <- lapply(x, function(v) if (is.character(v)) factor(v) else v)
  # target of the first step of the model: y - v-hat of the linear Fay-Herriot fit
  lin <- .fit_fh_lm(object$df_est$direct[s], .design_matrix(f0, object$data)[s, , drop = FALSE],
                    object$df_est$vardir[s], object$method)
  r <- object$df_est$direct[s] - drop(.design_matrix(f0, object$data)[s, , drop = FALSE] %*% lin$beta)
  v <- lin$sigma2_v / (lin$sigma2_v + object$df_est$vardir[s]) * r
  target <- object$df_est$direct[s] - v

  bor <- .with_seed(seed, Boruta::Boruta(
    x = x, y = target, pValue = p_value, maxRuns = max_runs,
    getImp = Boruta::getImpRfZ, num.threads = 1
  ))
  decision <- data.frame(
    variable = names(bor$finalDecision),
    decision = as.character(bor$finalDecision),
    median_importance = Boruta::attStats(bor)$medianImp,
    row.names = NULL, stringsAsFactors = FALSE
  )
  keep <- decision$variable[decision$decision == "Confirmed"]
  if (length(keep) == 0) {
    cli::cli_warn("No covariate was confirmed; the model is returned unchanged.")
    keep <- vars
  }
  final_formula <- stats::as.formula(paste(deparse(f0[[2]]), "~", paste(keep, collapse = " + ")),
                                     env = environment(f0))
  model <- .refit_flexsae(object, final_formula)
  structure(
    list(model = model, formula = final_formula, selected = keep, method = "boruta",
         decision = decision),
    class = "flexsae_selection"
  )
}

#' @rdname select_flexsae
#' @param x A `"flexsae_selection"` object.
#' @param ... Not used.
#' @export
print.flexsae_selection <- function(x, ...) {
  cat("Covariate selection (", x$method, ")\n", sep = "")
  cat("Selected model: ", deparse(x$formula), "\n", sep = "")
  if (x$method == "step") {
    cat("Criterion: ", x$criterion, ", direction: ", x$direction, "\n\n", sep = "")
    print(x$path, row.names = FALSE)
  } else {
    cat("\n")
    print(x$decision, row.names = FALSE)
  }
  invisible(x)
}
