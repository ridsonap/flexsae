#' Benchmark Area Estimates to Aggregate Targets
#'
#' Adjusts the area estimates of a fitted model so that their weighted
#' aggregate, overall or within groups, equals a target, and estimates the MSE
#' of the benchmarked estimates by applying the same benchmarking to every
#' bootstrap replicate of the model.
#'
#' @param object A `"flexsae"` model, or for the default method a numeric
#'   vector of estimates.
#' @param target Benchmark target(s) \eqn{T_g}: a number (the same for every
#'   group), a named vector by group, a data frame with columns `group` and
#'   `target`, or `"direct"` for the weighted aggregate of the direct estimates
#'   of the same survey. May be `NULL` when `national_target` is given (see
#'   Details).
#' @param weight Area weights \eqn{w_d} (e.g. population sizes): a column name of
#'   the model data, a one-sided formula or a vector. `NULL` gives equal
#'   weights. Must be positive.
#' @param group Optional groups (e.g. provinces), given like `weight`. Each
#'   group is benchmarked to its own target.
#' @param method Benchmarking rule (see Details): `"ratio"` (default),
#'   `"difference"`, `"optimal"` or `"logit"`.
#' @param type `"mean"` (default): the target is the weighted mean
#'   \eqn{\sum_d w_d\theta_d / \sum_d w_d} of the group; `"total"`: the target is
#'   \eqn{\sum_d w_d\theta_d}.
#' @param target_var Optional sampling variance of a numeric `target` (a
#'   number or one per group), when the target is itself an estimate from
#'   another source, independent of the direct estimates. `NULL` (default)
#'   treats a numeric target as known exactly.
#' @param national_target Optional overall target for two-stage benchmarking:
#'   the group targets (or, when `target = NULL`, the group aggregates of the
#'   estimates) are first adjusted with the same rule so that they aggregate to
#'   `national_target`; then the areas are benchmarked to them. Treated as
#'   known exactly.
#' @param areas `"all"` (default): every area enters the aggregate and is
#'   adjusted. `"sampled"`: only sampled areas; unsampled areas are left as
#'   they are.
#' @param bounds Optional `c(lower, upper)`; benchmarked estimates outside are
#'   flagged in column `out_of_bounds`, with a warning. They are not changed.
#' @param mse `"bootstrap"` (default) or `"none"`. The bootstrap needs a model
#'   fitted with `mse = "bootstrap"`.
#' @param seed Optional seed for the target errors drawn when `target_var` is
#'   given. The random number stream of the session is restored afterwards.
#' @param area_mse For the default method: optional MSE of the estimates, used
#'   by `method = "optimal"`.
#' @param ... Not used.
#'
#' @details
#' **Rules.** Within group \eqn{g}, with \eqn{a_d = w_d/\sum_{j \in g} w_j}
#' (`type = "mean"`) or \eqn{a_d = w_d} (`"total"`) and
#' \eqn{\delta_g = T_g - \sum_{d \in g} a_d\hat\theta_d}, every rule except
#' `"logit"` has the form \eqn{\tilde\theta_d = \hat\theta_d + k_d\delta_g} with
#' \eqn{\sum_d a_d k_d = 1}:
#' * `"ratio"`: \eqn{k_d = \hat\theta_d/\sum_j a_j\hat\theta_j}, i.e.
#'   \eqn{\tilde\theta_d = \hat\theta_d T_g/\sum_j a_j\hat\theta_j};
#' * `"difference"`: \eqn{k_d = 1/\sum_j a_j};
#' * `"optimal"`: \eqn{k_d = a_d\,\mathrm{mse}_d/\sum_j a_j^2\mathrm{mse}_j}, which
#'   minimizes \eqn{\sum_d(\tilde\theta_d - \hat\theta_d)^2/\mathrm{mse}_d} under the
#'   constraint (Datta et al., 2011; Bell, Datta and Ghosh, 2013);
#' * `"logit"`: \eqn{\tilde\theta_d = \mathrm{logit}^{-1}(\mathrm{logit}\,\hat\theta_d + c_g)}
#'   with \eqn{c_g} solving the constraint, for estimates in (0, 1).
#'
#' The rules agree with `benchmark_sae()` of \pkg{fastsae}; `"ratio"` and
#' `"difference"` agree with `type = "ratio"` and `"raking"` of
#' `emdi::benchmark()`, and `"optimal"` with `method = "optimum"` of
#' `saebenchmarking::sae_benchmarking()`. `emdi::benchmark(type = "MSE_adj")`
#' uses \eqn{k_d \propto \mathrm{mse}_d} instead of \eqn{a_d\,\mathrm{mse}_d}.
#' Two-stage benchmarking follows `fastsae::benchmark_sae()`; see Ghosh and
#' Steorts (2013) and Pfeffermann, Sikov and Tiller (2014) for the setting.
#'
#' **MSE.** The model bootstrap stores, for each replicate \eqn{b}, the
#' bootstrap truth \eqn{\theta^*}, the direct estimates \eqn{y^*} and the
#' refitted estimates \eqn{\hat\theta^*}. Each replicate is benchmarked with the
#' same rule to its bootstrap target, and
#' \eqn{\mathrm{mse}^{bm}_d = B^{-1}\sum_b(\tilde\theta^*_d - \theta^*_d)^2}, paired
#' with \eqn{\mathrm{mse}_d = B^{-1}\sum_b(\hat\theta^*_d - \theta^*_d)^2} over the
#' same replicates. The bootstrap target is
#' * known target: \eqn{T^*_g = \sum_d a_d\theta^*_d};
#' * `target_var` given: \eqn{T^*_g = \sum_d a_d\theta^*_d + \varepsilon^*_g},
#'   \eqn{\varepsilon^*_g \sim N(0, \mathrm{target\_var}_g)};
#' * `target = "direct"`: \eqn{T^*_g = \sum_d a_d y^*_d}, which keeps the
#'   correlation between the target and the estimates (the scheme of
#'   `saebenchmarking::mse_benchmarking()`, after Steorts and Ghosh, 2013).
#'
#' Because the model is refitted in every replicate, the MSE includes the
#' estimation of the model; `"optimal"` keeps \eqn{\mathrm{mse}_d} of the original
#' fit in all replicates. Benchmarking can raise or lower the MSE of an area;
#' `mse_ratio` is `mse_benchmarked / mse`. The bootstrap is done under the
#' fitted model, so it does not show how benchmarking reduces a bias of a
#' misspecified model.
#'
#' @return A data frame of class `c("flexsae_benchmark", "data.frame")` with
#'   columns `domain`, `group`, `weight`, `sampled`, `original`, `benchmarked`,
#'   `adjustment`, `rel_adjustment_pct`, `target`, `mse`, `mse_benchmarked`,
#'   `mse_ratio` and `out_of_bounds`, and attributes `method`, `type`,
#'   `target_source` (`"known"`, `"estimated"` or `"direct"`), `areas`,
#'   `hierarchical`, `verification` (aggregates by group), `national` (overall
#'   aggregate, two-stage only), `mse_method`, `B` and `n_failed`.
#'
#' @references
#' Bell, W. R., Datta, G. S. and Ghosh, M. (2013). Benchmarking small area
#' estimators. *Biometrika*, 100(1), 189-202. \doi{10.1093/biomet/ass063}
#'
#' Datta, G. S., Ghosh, M., Steorts, R. and Maples, J. (2011). Bayesian
#' benchmarking with applications to small area estimation. *TEST*, 20(3),
#' 574-588. \doi{10.1007/s11749-010-0218-y}
#'
#' Ghosh, M. and Steorts, R. C. (2013). Two-stage benchmarking as applied to
#' small area estimation. *TEST*, 22(4), 670-687.
#' \doi{10.1007/s11749-013-0338-2}
#'
#' Pfeffermann, D., Sikov, A. and Tiller, R. (2014). Single- and two-stage
#' cross-sectional and time series benchmarking procedures for small area
#' estimation. *TEST*, 23(4), 631-666. \doi{10.1007/s11749-014-0398-y}
#'
#' Steorts, R. and Ghosh, M. (2013). On estimation of mean squared errors of
#' benchmarked empirical Bayes estimators. *Statistica Sinica*, 23(2),
#' 749-767. \doi{10.5705/ss.2012.053}
#'
#' @export
#' @examples
#' set.seed(1)
#' D <- 40
#' dat <- data.frame(area = paste0("A", 1:D), x1 = runif(D, -2, 2),
#'                   prov = rep(c("P1", "P2"), each = D / 2), N = rpois(D, 500))
#' dat$vardir <- runif(D, 0.2, 0.6)
#' dat$y <- 5 + sin(1.5 * dat$x1) + rnorm(D, sd = 0.5) + rnorm(D, sd = sqrt(dat$vardir))
#'
#' fit <- pspline_sae(y ~ s(x1, k = 5), vardir = "vardir", domain = "area",
#'                    data = dat, B = 20, seed = 1)
#'
#' # one target per province, known exactly
#' bm <- benchmark_flexsae(fit, target = c(P1 = 5.1, P2 = 4.9), weight = "N",
#'                         group = "prov")
#' bm
#'
#' # self-consistency with the direct estimates of the same survey
#' bm_dir <- benchmark_flexsae(fit, target = "direct", weight = "N",
#'                             method = "difference")
#' summary(bm_dir)
benchmark_flexsae <- function(object, ...) {
  UseMethod("benchmark_flexsae")
}

#' @rdname benchmark_flexsae
#' @export
benchmark_flexsae.flexsae <- function(
  object,
  target = NULL,
  weight = NULL,
  group = NULL,
  method = c("ratio", "difference", "optimal", "logit"),
  type = c("mean", "total"),
  target_var = NULL,
  national_target = NULL,
  areas = c("all", "sampled"),
  bounds = NULL,
  mse = c("bootstrap", "none"),
  seed = NULL,
  ...
) {
  method <- match.arg(method)
  type <- match.arg(type)
  areas <- match.arg(areas)
  mse <- match.arg(mse)
  .check_no_dots(...)
  df <- object$df_est
  D <- nrow(df)
  w <- if (is.null(weight)) rep(1, D) else .get_variable(object$data, weight, "weight")
  grp <- if (is.null(group)) rep("All", D) else .get_variable(object$data, group, "group")

  boot <- NULL
  if (mse == "bootstrap") {
    if (is.null(object$boot_theta) || is.null(object$boot_y)) {
      cli::cli_abort(c(
        "The model has no bootstrap replicates.",
        "i" = "Refit it with {.code mse = \"bootstrap\"}, or use {.code mse = \"none\"}."
      ), call = NULL)
    }
    boot <- list(est = object$boot_errors + object$boot_theta, theta = object$boot_theta, y = object$boot_y)
  }
  .benchmark_core(
    domain = df$domain, est = df$est, direct = df$direct, sampled = df$sampled,
    area_mse = df$mse, weight = w, group = grp, target = target, target_var = target_var,
    national_target = national_target, method = method, type = type, areas = areas,
    bounds = bounds, boot = boot, seed = seed
  )
}

#' @rdname benchmark_flexsae
#' @export
benchmark_flexsae.default <- function(
  object,
  target = NULL,
  weight = NULL,
  group = NULL,
  method = c("ratio", "difference", "optimal", "logit"),
  type = c("mean", "total"),
  target_var = NULL,
  national_target = NULL,
  bounds = NULL,
  area_mse = NULL,
  ...
) {
  method <- match.arg(method)
  type <- match.arg(type)
  .check_no_dots(...)
  if (!is.numeric(object) || anyNA(object)) {
    cli::cli_abort("{.arg object} must be a numeric vector without missing values or a {.cls flexsae} model.", call = NULL)
  }
  D <- length(object)
  if (identical(target, "direct")) {
    cli::cli_abort("{.code target = \"direct\"} needs a {.cls flexsae} model.", call = NULL)
  }
  df <- data.frame(row = seq_len(D))
  w <- if (is.null(weight)) rep(1, D) else .get_variable(df, weight, "weight")
  grp <- if (is.null(group)) rep("All", D) else .get_variable(df, group, "group")
  if (!is.null(area_mse) && (!is.numeric(area_mse) || length(area_mse) != D)) {
    cli::cli_abort("{.arg area_mse} must be a numeric vector of length {D}.", call = NULL)
  }
  .benchmark_core(
    domain = names(object) %||% seq_len(D), est = as.numeric(object), direct = rep(NA_real_, D),
    sampled = rep(TRUE, D), area_mse = area_mse %||% rep(NA_real_, D), weight = w,
    group = grp, target = target, target_var = target_var, national_target = national_target,
    method = method, type = type, areas = "all", bounds = bounds, boot = NULL, seed = NULL
  )
}

#' Shared benchmarking worker
#' @noRd
.benchmark_core <- function(domain, est, direct, sampled, area_mse, weight, group, target,
                            target_var, national_target, method, type, areas, bounds, boot, seed) {
  D <- length(est)
  if (!is.numeric(weight) || anyNA(weight) || any(!is.finite(weight)) || any(weight <= 0)) {
    cli::cli_abort("{.arg weight} must be positive and finite for every area.", call = NULL)
  }
  if (anyNA(group)) {
    cli::cli_abort("{.arg group} must not contain missing values.", call = NULL)
  }
  group <- as.character(group)
  if (!is.null(national_target) &&
    (!is.numeric(national_target) || length(national_target) != 1 || !is.finite(national_target))) {
    cli::cli_abort("{.arg national_target} must be a single finite number.", call = NULL)
  }
  if (!is.null(bounds) && (!is.numeric(bounds) || length(bounds) != 2 || anyNA(bounds) || bounds[1] >= bounds[2])) {
    cli::cli_abort("{.arg bounds} must be {.code c(lower, upper)} with lower < upper.", call = NULL)
  }

  agg <- if (areas == "all") rep(TRUE, D) else sampled
  groups <- unique(group[agg])
  hierarchical <- !is.null(national_target) && length(groups) > 1
  wa <- weight[agg]
  ga <- group[agg]
  a <- .agg_weights(wa, ga, type)

  # targets
  direct_target <- identical(target, "direct")
  if (direct_target) {
    if (!is.null(target_var)) {
      cli::cli_abort("{.arg target_var} cannot be used with {.code target = \"direct\"}.", call = NULL)
    }
    if (any(!sampled[agg])) {
      cli::cli_abort(c(
        "{.code target = \"direct\"} needs a direct estimate for every area in the aggregate.",
        "i" = "Use {.code areas = \"sampled\"} to benchmark the sampled areas only."
      ), call = NULL)
    }
    target_g <- tapply(a * direct[agg], factor(ga, levels = groups), sum)
  } else {
    target_g <- .resolve_targets(target, groups, national_target, hierarchical)
  }
  target_var_g <- .resolve_targets(target_var, groups, NULL, FALSE, arg = "target_var", allow_null = TRUE)
  if (!is.null(target_var_g)) {
    if (is.null(target)) {
      cli::cli_abort("{.arg target_var} needs a numeric {.arg target}.", call = NULL)
    }
    if (any(target_var_g < 0)) cli::cli_abort("{.arg target_var} must not be negative.", call = NULL)
  }
  target_source <- if (direct_target) "direct" else if (is.null(target_var_g)) "known" else "estimated"

  mse_a <- area_mse[agg]
  if (method == "optimal" && (anyNA(mse_a) || any(mse_a < 0))) {
    cli::cli_abort(c(
      "{.code method = \"optimal\"} needs the MSE of every area in the aggregate.",
      "i" = "Fit the model with {.code mse = \"bootstrap\"}."
    ), call = NULL)
  }
  group_mse <- NULL
  if (method == "optimal" && hierarchical) {
    group_mse <- if (!is.null(boot)) {
      e <- boot$est[agg, , drop = FALSE] - boot$theta[agg, , drop = FALSE]
      rowMeans(rowsum(a * e, factor(ga, levels = groups))^2)
    } else {
      as.numeric(tapply(a^2 * mse_a, factor(ga, levels = groups), sum))
    }
    names(group_mse) <- groups
  }

  rule <- function(e, tg) {
    .bm_rule(e, wa, a, ga, groups, tg, national_target, method, type, mse_a, group_mse)
  }
  res <- rule(est[agg], target_g)
  benchmarked <- est
  benchmarked[agg] <- res$bench

  # bootstrap MSE: benchmark every replicate to its bootstrap target
  mse_orig <- area_mse
  mse_bm <- rep(NA_real_, D)
  n_failed <- NA_integer_
  B <- NA_integer_
  if (!is.null(boot)) {
    B <- ncol(boot$theta)
    err <- boot$est - boot$theta
    err_bm <- err
    ok <- rep(TRUE, B)
    fac <- factor(ga, levels = groups)
    .with_seed(seed, {
      for (b in seq_len(B)) {
        th <- boot$theta[agg, b]
        tg_star <- if (is.null(target) && hierarchical) {
          NULL
        } else if (direct_target) {
          as.numeric(tapply(a * boot$y[agg, b], fac, sum))
        } else {
          tt <- as.numeric(tapply(a * th, fac, sum))
          if (!is.null(target_var_g)) tt <- tt + stats::rnorm(length(groups), 0, sqrt(target_var_g))
          tt
        }
        if (!is.null(tg_star)) names(tg_star) <- groups
        nat_star <- if (hierarchical) .national_aggregate(th, wa, type) else NULL
        bb <- tryCatch(
          .bm_rule(boot$est[agg, b], wa, a, ga, groups, tg_star, nat_star, method, type, mse_a, group_mse)$bench,
          error = function(e) NULL
        )
        if (is.null(bb) || anyNA(bb)) {
          ok[b] <- FALSE
          next
        }
        err_bm[agg, b] <- bb - th
      }
    })
    n_failed <- sum(!ok)
    if (n_failed == B) {
      cli::cli_abort("The benchmarking failed in every bootstrap replicate.", call = NULL)
    }
    if (n_failed > 0) {
      cli::cli_warn("Benchmarking failed in {n_failed} of {B} bootstrap replicate{?s}; the MSE uses the others.")
    }
    mse_orig <- rowMeans(err[, ok, drop = FALSE]^2)
    mse_bm <- rowMeans(err_bm[, ok, drop = FALSE]^2)
  }

  out_of_bounds <- if (is.null(bounds)) rep(NA, D) else benchmarked < bounds[1] | benchmarked > bounds[2]
  if (isTRUE(any(out_of_bounds))) {
    cli::cli_warn("{sum(out_of_bounds)} benchmarked estimate{?s} {?is/are} outside {.arg bounds}.")
  }
  target_col <- rep(NA_real_, D)
  target_col[agg] <- res$target[ga]

  out <- data.frame(
    domain = domain,
    group = group,
    weight = weight,
    sampled = sampled,
    original = est,
    benchmarked = benchmarked,
    adjustment = benchmarked - est,
    rel_adjustment_pct = ifelse(abs(est) > .Machine$double.eps, 100 * (benchmarked - est) / est, NA_real_),
    target = target_col,
    mse = if (is.null(boot)) area_mse else mse_orig,
    mse_benchmarked = mse_bm,
    mse_ratio = mse_bm / (if (is.null(boot)) area_mse else mse_orig),
    out_of_bounds = out_of_bounds,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  fac <- factor(ga, levels = groups)
  verification <- data.frame(
    group = groups,
    target = unname(res$target[groups]),
    original = as.numeric(tapply(a * est[agg], fac, sum)),
    benchmarked = as.numeric(tapply(a * res$bench, fac, sum)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  verification$discrepancy <- verification$benchmarked - verification$target
  national <- NULL
  if (hierarchical) {
    national <- data.frame(
      target = national_target,
      original = .national_aggregate(est[agg], wa, type),
      benchmarked = .national_aggregate(res$bench, wa, type)
    )
    national$discrepancy <- national$benchmarked - national$target
  }

  structure(
    out,
    class = c("flexsae_benchmark", "data.frame"),
    method = method,
    type = type,
    target_source = target_source,
    areas = areas,
    hierarchical = hierarchical,
    verification = verification,
    national = national,
    mse_method = if (is.null(boot)) "none" else "bootstrap",
    B = B,
    n_failed = n_failed
  )
}

#' Aggregation weights: normalized within group for type "mean"
#' @noRd
.agg_weights <- function(w, g, type) {
  if (type == "total") {
    return(w)
  }
  w / stats::ave(w, g, FUN = sum)
}

#' Overall aggregate across groups (two-stage benchmarking)
#' @noRd
.national_aggregate <- function(x, w, type) {
  if (type == "mean") sum(w * x) / sum(w) else sum(w * x)
}

#' Apply a benchmarking rule (stage 1 for two-stage, then within groups)
#'
#' Follows `.benchmark_worker()` of fastsae.
#' @param est,w,a,g Estimates, weights, aggregation weights and groups of the
#'   aggregated areas.
#' @param target_g Named group targets, or NULL (two-stage: group aggregates).
#' @return List with `bench` and the final group targets `target`.
#' @noRd
.bm_rule <- function(est, w, a, g, groups, target_g, national, method, type, mse, group_mse) {
  fac <- factor(g, levels = groups)
  agg_g <- tapply(a * est, fac, sum)
  if (is.null(target_g)) target_g <- agg_g
  target_g <- target_g[groups]

  if (!is.null(national) && length(groups) == 1) {
    target_g[] <- national
  } else if (!is.null(national)) {
    # stage 1: group targets to the national target, with group shares W_g
    W <- as.numeric(tapply(w, fac, sum)) / sum(w)
    nat_now <- if (type == "mean") sum(W * target_g) else sum(target_g)
    diff_nat <- national - nat_now
    if (abs(diff_nat) > 1e-7) {
      n_g <- length(groups)
      target_g <- switch(method,
        ratio = target_g * national / nat_now,
        difference = target_g + if (type == "mean") diff_nat else diff_nat / n_g,
        optimal = {
          m_g <- pmax(group_mse[groups], 1e-8)
          a_g <- if (type == "mean") W else rep(1, n_g)
          target_g + a_g * m_g * diff_nat / sum(a_g^2 * m_g)
        },
        logit = {
          if (type != "mean") {
            cli::cli_abort("Two-stage {.val logit} benchmarking needs {.code type = \"mean\"}.", call = NULL)
          }
          lt <- stats::qlogis(target_g)
          root <- stats::uniroot(function(c) sum(W * stats::plogis(lt + c)) - national,
                                 c(-50, 50), tol = 1e-12)$root
          stats::plogis(lt + root)
        }
      )
      names(target_g) <- groups
    }
  }

  bench <- est
  for (gr in groups) {
    i <- which(g == gr)
    ai <- a[i]
    ei <- est[i]
    tg <- target_g[[gr]]
    agg <- sum(ai * ei)
    bench[i] <- switch(method,
      ratio = {
        if (abs(agg) < .Machine$double.eps) {
          cli::cli_abort("The aggregate of group {.val {gr}} is zero; ratio benchmarking is undefined.", call = NULL)
        }
        ei * tg / agg
      },
      difference = ei + (tg - agg) / sum(ai),
      optimal = {
        m <- pmax(mse[i], 1e-8)
        ei + ai * m * (tg - agg) / sum(ai^2 * m)
      },
      logit = {
        if (any(ei <= 0 | ei >= 1)) {
          cli::cli_abort("Logit benchmarking needs estimates in (0, 1); group {.val {gr}} has some outside.", call = NULL)
        }
        if (type == "mean" && (tg <= 0 || tg >= 1)) {
          cli::cli_abort("Logit benchmarking with {.code type = \"mean\"} needs targets in (0, 1).", call = NULL)
        }
        le <- stats::qlogis(ei)
        f <- function(c) sum(ai * stats::plogis(le + c)) - tg
        lim <- if (f(-20) * f(20) > 0) 50 else 20
        stats::plogis(le + stats::uniroot(f, c(-lim, lim), tol = 1e-12)$root)
      }
    )
  }
  list(bench = bench, target = target_g)
}

#' Resolve targets given as a number, named vector or data frame
#' @noRd
.resolve_targets <- function(target, groups, national_target, hierarchical,
                             arg = "target", allow_null = FALSE) {
  if (is.null(target)) {
    if (allow_null || hierarchical) {
      return(NULL)
    }
    if (!is.null(national_target)) {
      return(stats::setNames(rep(national_target, length(groups)), groups))
    }
    cli::cli_abort("Give {.arg target} or {.arg national_target}.", call = NULL)
  }
  if (is.data.frame(target)) {
    if (!all(c("group", "target") %in% names(target))) {
      cli::cli_abort("A data frame {.arg {arg}} needs columns {.val group} and {.val target}.", call = NULL)
    }
    target <- stats::setNames(target$target, as.character(target$group))
  }
  if (!is.numeric(target)) {
    cli::cli_abort("{.arg {arg}} must be numeric, a named vector, a data frame or {.val direct}.", call = NULL)
  }
  out <- if (length(target) == 1 && is.null(names(target))) {
    rep(as.numeric(target), length(groups))
  } else if (!is.null(names(target))) {
    miss <- setdiff(groups, names(target))
    if (length(miss) > 0) {
      cli::cli_abort("{.arg {arg}} has no value for group{?s} {.val {miss}}.", call = NULL)
    }
    as.numeric(target[groups])
  } else if (length(target) == length(groups)) {
    as.numeric(target)
  } else {
    cli::cli_abort("{.arg {arg}} must have one value per group ({length(groups)}), not {length(target)}.", call = NULL)
  }
  if (any(!is.finite(out))) {
    cli::cli_abort("{.arg {arg}} must be finite.", call = NULL)
  }
  stats::setNames(out, groups)
}

#' Abort when unused arguments are passed through ...
#' @noRd
.check_no_dots <- function(...) {
  if (...length() > 0) {
    nms <- names(list(...)) %||% rep("", ...length())
    nms[!nzchar(nms)] <- "(unnamed)"
    cli::cli_abort("Unknown argument{?s}: {.val {nms}}.", call = NULL)
  }
  invisible(TRUE)
}

#' @rdname benchmark_flexsae
#' @param x A `"flexsae_benchmark"` object.
#' @param n Number of areas to print.
#' @export
print.flexsae_benchmark <- function(x, n = 6, ...) {
  # `[` keeps the class but drops the attributes: print a subset as a data frame
  if (is.null(attr(x, "method"))) {
    print(as.data.frame(x), ...)
    return(invisible(x))
  }
  src <- switch(attr(x, "target_source"),
    known = "known exactly",
    estimated = "estimated (with target_var)",
    direct = "direct estimates of the same survey"
  )
  cat("Benchmarked area estimates\n")
  cat("Method: ", attr(x, "method"), ", type: ", attr(x, "type"), ", areas: ", attr(x, "areas"), "\n", sep = "")
  cat("Target: ", src, if (isTRUE(attr(x, "hierarchical"))) ", two-stage" else "", "\n", sep = "")
  nat <- attr(x, "national")
  if (!is.null(nat)) {
    cat(sprintf("National: target %.6g, benchmarked aggregate %.6g\n", nat$target, nat$benchmarked))
  }
  v <- attr(x, "verification")
  cat("\nAggregates by group:\n")
  print(utils::head(v, 10), row.names = FALSE, digits = 6)
  if (nrow(v) > 10) cat("... and", nrow(v) - 10, "more groups\n")
  if (identical(attr(x, "mse_method"), "bootstrap")) {
    cat(sprintf("\nMSE: bootstrap, %d replicates (%d failed)\n", attr(x, "B") - attr(x, "n_failed"), attr(x, "n_failed")))
  }
  cat("\n")
  print(utils::head(as.data.frame(x), n), row.names = FALSE, digits = 4)
  if (nrow(x) > n) cat("... and", nrow(x) - n, "more areas\n")
  invisible(x)
}

#' @rdname benchmark_flexsae
#' @export
summary.flexsae_benchmark <- function(object, ...) {
  print(object, n = 0)
  cat("Adjustment:\n")
  print(summary(object$adjustment))
  cat("\nRelative adjustment (%):\n")
  print(summary(object$rel_adjustment_pct))
  if (!all(is.na(object$mse_benchmarked))) {
    r <- object$mse_ratio
    cat("\nMSE benchmarked / MSE:\n")
    print(summary(r))
    cat(sprintf("MSE higher after benchmarking in %d of %d areas\n", sum(r > 1, na.rm = TRUE), sum(!is.na(r))))
  }
  if (!all(is.na(object$out_of_bounds))) {
    cat(sprintf("Outside bounds: %d\n", sum(object$out_of_bounds)))
  }
  invisible(object)
}
