#' Covariate Selection for flexsae Models
#'
#' Selects the covariates of a fitted model and refits it with them, or, for
#' [scfh_sae()], the number of clusters and the spatial penalty.
#'
#' @param object A model from [pspline_sae()], [merf_sae()] or [scfh_sae()].
#' @param method `"step"` (stepwise search on an information criterion;
#'   `pspline_sae()` models), `"boruta"` (all-relevant selection with the
#'   Boruta algorithm; `merf_sae()` models) or `"grid"` (\eqn{K} and
#'   \eqn{\phi}; `scfh_sae()` models). The default depends on the model.
#' @param criterion For `"step"` and `"grid"`: `"BIC"` (default), `"AIC"` or
#'   `"KIC"`.
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
#' @param K,phi For `"grid"`: candidate numbers of clusters and spatial
#'   penalties.
#' @param band For `"grid"`: range of the moderate penalties used by the
#'   selection rule.
#' @param delta For `"grid"`: tolerance of the stability plateau.
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
#' **Grid (`scfh_sae()`).** The two-step rule of Maranzano, Mattera and
#' Sugasawa (2026, sec. 4.3). Every pair of `K` and `phi` is fitted (without
#' MSE; `K = 1`, the Fay-Herriot model, once); a fit is admissible when every
#' cluster has enough sampled areas. The criterion uses the unpenalized
#' log-likelihood \eqn{\ell} of the sampled areas at the fitted partition and
#' \eqn{k = K(P + 2)} parameters (\eqn{P + 1} coefficients and one variance
#' per cluster), with the formulas above. (1) For every penalty in `band`, the
#' admissible `K` with the smallest criterion; \eqn{\hat K} is the most
#' frequent one, ties to the smaller. (2) At \eqn{\hat K}, the stability
#' \eqn{S(\phi)} is the mean adjusted Rand index (Hubert and Arabie, 1985)
#' between the partition at \eqn{\phi} and those at the neighbouring
#' admissible penalties of the grid in `band` (our reading of "adjacent grid
#' values"); \eqn{\hat\phi} is the smallest penalty
#' with \eqn{S(\phi) \ge \max S - \delta}. The paper gives the number of
#' parameters as \eqn{K(P + 2)} without a fuller definition of its "effective
#' degrees of freedom". The partition is optimized on the data, so in-sample
#' criteria still tend to choose too many clusters (sec. 4.3 of the paper);
#' check the grid and the Fay-Herriot model (`K = 1`) as well. The selected
#' model is refitted with the MSE settings of `object`.
#'
#' @return An object of class `"flexsae_selection"`: a list with `model` (the
#'   refitted model), `formula`, `selected` (kept terms or variables; for
#'   `"grid"` the selected `K` and `phi`), `method`, and for `"step"` the
#'   `path` (data frame of steps) and `criterion`, for `"boruta"` the
#'   `decision` table, for `"grid"` the `grid` table (one row per fit, with
#'   `stability` at \eqn{\hat K}) and `criterion`.
#'
#' @references
#' Hubert, L. and Arabie, P. (1985). Comparing partitions. *Journal of
#' Classification*, 2(1), 193-218. \doi{10.1007/BF01908075}
#'
#' Kursa, M. B. and Rudnicki, W. R. (2010). Feature selection with the Boruta
#' package. *Journal of Statistical Software*, 36(11), 1-13.
#' \doi{10.18637/jss.v036.i11}
#'
#' Marhuenda, Y., Morales, D. and Pardo, M. C. (2014). Information criteria
#' for Fay-Herriot model selection. *Computational Statistics & Data
#' Analysis*, 70, 268-280. \doi{10.1016/j.csda.2013.09.016}
#'
#' Maranzano, P., Mattera, R. and Sugasawa, S. (2026). Small area estimation
#' under spatial regimes: spatially clustered Fay-Herriot models for
#' agricultural indicators. arXiv preprint 2608.13638.
#' \doi{10.48550/arXiv.2608.13638}
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
  seed = NULL,
  K = 1:4,
  phi = seq(0, 1, by = 0.125),
  band = c(0.25, 1),
  delta = 0.05
) {
  if (!inherits(object, "flexsae")) {
    cli::cli_abort("{.arg object} must be a model from {.fn pspline_sae}, {.fn merf_sae} or {.fn scfh_sae}.",
                   call = NULL)
  }
  if (identical(object$model, "boost_sae")) {
    cli::cli_abort(c(
      "Covariate selection for {.fn boost_sae} models is not implemented.",
      "i" = "See {.fn importance_flexsae} for the gain of each covariate."
    ), call = NULL)
  }
  is_merf <- identical(object$model, "merf_sae")
  is_scfh <- identical(object$model, "scfh_sae")
  method <- method %||% if (is_merf) "boruta" else if (is_scfh) "grid" else "step"
  method <- match.arg(method, c("step", "boruta", "grid"))
  if ((method == "grid") != is_scfh) {
    cli::cli_abort(c(
      if (is_scfh) "{.fn scfh_sae} models are tuned with {.code method = \"grid\"}." else
        "{.code method = \"grid\"} is available for {.fn scfh_sae} models.",
      "i" = "Covariate selection for {.fn scfh_sae} is not implemented."
    ), call = NULL)
  }
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
  if (method == "grid") {
    .select_grid(object, criterion, K, phi, band, delta, trace)
  } else if (method == "step") {
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

#' Grid search over K and phi for scfh_sae()
#' @noRd
.select_grid <- function(object, criterion, K, phi, band, delta, trace) {
  K <- sort(unique(K))
  phi <- sort(unique(phi))
  if (!is.numeric(K) || length(K) == 0 || any(!is.finite(K) | K < 1 | K != round(K))) {
    cli::cli_abort("{.arg K} must be positive whole numbers.", call = NULL)
  }
  if (!is.numeric(phi) || length(phi) == 0 || any(!is.finite(phi) | phi < 0)) {
    cli::cli_abort("{.arg phi} must be non-negative numbers.", call = NULL)
  }
  if (!is.numeric(band) || length(band) != 2 || band[1] > band[2]) {
    cli::cli_abort("{.arg band} must be two increasing numbers.", call = NULL)
  }
  in_band <- phi[phi >= band[1] & phi <= band[2]]
  if (length(in_band) == 0) {
    cli::cli_abort("No value of {.arg phi} lies in {.arg band} [{band[1]}, {band[2]}].", call = NULL)
  }
  if (!is.numeric(delta) || length(delta) != 1 || !is.finite(delta) || delta < 0) {
    cli::cli_abort("{.arg delta} must be a single non-negative number.", call = NULL)
  }
  n_coef <- ncol(.design_matrix(object$formula, object$data))
  m <- sum(object$df_est$sampled)
  ic <- function(ll, k) switch(criterion, AIC = -2 * ll + 2 * k, BIC = -2 * ll + k * log(m), KIC = -2 * ll + 3 * k)

  grid <- expand.grid(phi = phi, K = K)[, c("K", "phi")]
  grid <- grid[!(grid$K == 1 & grid$phi != phi[1]), ]
  fits <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    fits[[i]] <- tryCatch(
      suppressWarnings(.refit_scfh(object, grid$K[i], grid$phi[i], mse = "none")),
      error = function(e) NULL
    )
  }
  ok <- !vapply(fits, is.null, logical(1))
  grid$loglik <- vapply(fits, function(f) if (is.null(f)) NA_real_ else f$loglik, numeric(1))
  grid$admissible <- ok & vapply(fits, function(f) isTRUE(f$admissible), logical(1))
  grid$criterion <- ic(grid$loglik, grid$K * (n_coef + 1))
  grid$criterion[!grid$admissible] <- NA_real_
  if (trace) print(grid, row.names = FALSE)

  # (1) K: most frequent minimizer over the band; K = 1 is in every comparison
  pick <- vapply(in_band, function(p) {
    rows <- which((grid$phi == p | grid$K == 1) & grid$admissible)
    if (length(rows) == 0) NA_integer_ else as.integer(grid$K[rows][which.min(grid$criterion[rows])])
  }, integer(1))
  if (all(is.na(pick))) {
    cli::cli_abort("No admissible fit in {.arg band}; use smaller values of {.arg K}.", call = NULL)
  }
  freq <- table(pick)
  K_hat <- min(as.integer(names(freq)[freq == max(freq)]))

  # (2) phi: smallest penalty on the stability plateau at K_hat
  grid$stability <- NA_real_
  if (K_hat == 1) {
    phi_hat <- phi[1]
  } else {
    rows <- which(grid$K == K_hat & grid$phi %in% in_band & grid$admissible)
    s <- object$df_est$sampled
    for (j in seq_along(rows)) {
      nbr <- rows[abs(j - seq_along(rows)) == 1]
      grid$stability[rows[j]] <- if (length(nbr) == 0) 1 else mean(vapply(nbr, function(r) {
        .adjusted_rand(fits[[rows[j]]]$cluster$cluster[s], fits[[r]]$cluster$cluster[s])
      }, numeric(1)))
    }
    st <- grid$stability[rows]
    phi_hat <- min(grid$phi[rows][st >= max(st) - delta])
  }
  model <- .refit_scfh(object, K_hat, phi_hat, mse = object$mse_method)
  structure(
    list(model = model, formula = object$formula, selected = c(K = K_hat, phi = phi_hat),
         method = "grid", grid = grid, criterion = criterion, band = band),
    class = "flexsae_selection"
  )
}

#' Refit an scfh_sae() model with another K and phi
#' @noRd
.refit_scfh <- function(object, K, phi, mse = object$mse_method) {
  args <- list(
    formula = object$formula, vardir = object$df_est$vardir, data = object$data, W = object$W,
    K = K, phi = phi, method = object$method, n_start = object$n_start, maxit = object$maxit,
    mse = mse, seed = object$seed, conf_level = object$conf_level
  )
  if (mse == "bootstrap") args$B <- object$B
  fit <- do.call(scfh_sae, args)
  fit$df_est$domain <- fit$cluster$domain <- object$df_est$domain
  fit
}

#' @rdname select_flexsae
#' @param x A `"flexsae_selection"` object.
#' @param ... Not used.
#' @export
print.flexsae_selection <- function(x, ...) {
  if (x$method == "grid") {
    cat("Selection of K and phi (", x$criterion, ", stability plateau in [", x$band[1], ", ",
        x$band[2], "])\n", sep = "")
    cat("Selected: K = ", x$selected[["K"]], ", phi = ", format(x$selected[["phi"]]), "\n\n", sep = "")
    print(x$grid, row.names = FALSE)
    return(invisible(x))
  }
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
