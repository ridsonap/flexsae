#' Spatially Clustered Fay-Herriot Model
#'
#' @description
#' **Experimental.** Fits the
#' spatially clustered Fay-Herriot model (SC-FH) of Maranzano, Mattera and
#' Sugasawa (2026): the areas fall into \eqn{K} latent regimes with their own
#' regression coefficients and random effect variance,
#' \deqn{y_d = x_d'\theta_{k_d} + u_d + e_d,\quad u_d \sim N(0, \sigma^2_{u,k_d}),
#'   \quad e_d \sim N(0, \psi_d),}
#' and the partition is estimated jointly with the parameters by maximizing the
#' Fay-Herriot log-likelihood plus a Potts penalty
#' \eqn{\phi \sum_{d < d'} w_{dd'} 1(k_d = k_{d'})} that favours neighbouring
#' areas in the same regime.
#'
#' @param formula Two-sided formula of linear terms, e.g. `y ~ x1 + x2`.
#'   Rows with a missing response are unsampled areas.
#' @param vardir Sampling variances: column name, one-sided formula or vector.
#' @param domain Optional domain labels (column name, formula or vector).
#' @param data Data frame with one row per area.
#' @param W Square matrix over all rows of `data` (sampled and unsampled
#'   areas, in the same order), e.g. from `fastsae::create_weights()`. Two areas
#'   are neighbours when `W[d, d'] > 0` or `W[d', d] > 0`; the weights
#'   themselves are not used (a binary adjacency, as in the paper).
#' @param K Number of clusters. `K = 1` is the Fay-Herriot model.
#' @param phi Spatial penalty \eqn{\phi \ge 0}. `phi = 0` gives the
#'   non-spatial clusterwise Fay-Herriot model.
#' @param method Estimation of the cluster variances: `"AREML"` (default,
#'   adjusted REML of Li and Lahiri, 2010), `"REML"` or `"ML"`.
#' @param n_start Number of k-means restarts for the initial partition.
#' @param maxit Maximum number of alternating iterations.
#' @param mse `"bootstrap"` (default) or `"none"`.
#' @param B Number of bootstrap replicates.
#' @param seed Optional seed; the random number stream of the session is
#'   restored afterwards.
#' @param conf_level Confidence level of the intervals in `df_est`.
#'
#' @details
#' The algorithm follows Sections 2.3 and 4.1 of Maranzano et al. (2026):
#'
#' 1. Initial partition: k-means on the standardized covariates of the sampled
#'    areas, compressed by principal components (components up to 90% of the
#'    variance), keeping among `n_start` restarts the one with the smallest
#'    within-cluster sum of squares whose clusters all have at least
#'    \eqn{P + 2} areas (\eqn{P} covariates, so one more than the number of
#'    coefficients).
#' 2. Parameter step: a Fay-Herriot model per cluster. Adjusted REML maximizes
#'    \eqn{\log\sigma^2_u + \ell_R(\sigma^2_u)}, whose maximum is strictly
#'    positive, so that no cluster collapses to a synthetic predictor.
#' 3. Label step: one sequential sweep (iterated conditional modes, Besag,
#'    1986) over the sampled areas,
#'    \eqn{k_d \leftarrow \arg\max_k \{\log f(y_d; x_d'\theta_k, \sigma^2_{uk} + \psi_d) +
#'    \phi \sum_{d' \sim d} 1(k = k_{d'})\}}; the current label is kept on ties.
#'    A cluster with fewer than \eqn{P + 2} sampled areas keeps its last
#'    parameters.
#' 4. Stop when the partition does not change; if a partition repeats, the
#'    visited partition with the largest penalized objective is used. The
#'    clusters are then refitted once on the final partition and numbered by
#'    decreasing size.
#'
#' A covariate that is constant or collinear within a cluster (e.g. a factor
#' level absent from it) is dropped from that cluster's fit, as in [lm()]: its
#' coefficient is `NA` and it does not contribute to the cluster's predictions
#' (SCFH-6).
#'
#' The paper fits the sampled areas only. Here an unsampled area takes the
#' most frequent label among its labelled neighbours (the label step with the
#' penalty term only), propagated until no label changes; an area without
#' labelled neighbours takes the cluster whose mean of the standardized
#' covariates is closest. Unsampled areas get the synthetic estimate
#' \eqn{x_d'\hat\theta_{k_d}}. For a model without covariates the initial
#' k-means uses the direct estimates instead of the spatial coordinates of the
#' paper. These are our own choices (`technical-report/DISCREPANCIES.md`,
#' SCFH-1 and SCFH-2).
#'
#' The bootstrap MSE is Algorithm 2 of the paper: \eqn{\tau^*_d =
#' x_d'\hat\theta_{k_d} + u^*_d}, \eqn{u^*_d \sim N(0, \hat\sigma^2_{u,k_d})},
#' \eqn{y^*_d = \tau^*_d + e^*_d}, and the whole procedure, clustering
#' included, is run again on \eqn{y^*}, so the MSE carries the uncertainty of
#' the partition. The bootstrap labels are mapped to the clusters of the fit
#' by the majority rule and stored in `boot_cluster`.
#'
#' The label step is greedy and starts from a partition of the covariates, so
#' a large penalty can lock it into that initial partition: in our lattice
#' tests \eqn{\phi \ge 3} (rook neighbours) gave spatially smooth but wrong
#' regimes. Moderate values, as in the selection rule of the paper, are safer.
#'
#' The method comes from a preprint that had not been peer reviewed when it
#' was implemented, and no reference implementation was available; the
#' function is experimental. Choose `K` and `phi` with [select_flexsae()].
#'
#' @return An object of class `c("flexsae_scfh", "flexsae")` with `df_est`
#'   (standard table of estimates), `cluster` (data frame of domain and
#'   cluster), `coefficients` (by cluster, rows named `"[k] term"`),
#'   `sigma2_v` (one variance per cluster), `size`, `K`, `phi`, `loglik`
#'   (unpenalized log-likelihood of the sampled areas), `objective`
#'   (penalized), `iterations`, `admissible` (every cluster has at least
#'   \eqn{P + 2} sampled areas), `convergence`, the bootstrap replicates
#'   (`boot_errors`, `boot_theta`, `boot_y`, `boot_cluster`) and the settings.
#'
#' @references
#' Maranzano, P., Mattera, R. and Sugasawa, S. (2026). Small area estimation
#' under spatial regimes: spatially clustered Fay-Herriot models for
#' agricultural indicators. arXiv preprint 2608.13638.
#' \doi{10.48550/arXiv.2608.13638}
#'
#' Sugasawa, S. and Murakami, D. (2021). Spatially clustered regression.
#' *Spatial Statistics*, 44, 100525. \doi{10.1016/j.spasta.2021.100525}
#'
#' Li, H. and Lahiri, P. (2010). An adjusted maximum likelihood method for
#' solving small area estimation problems. *Journal of Multivariate Analysis*,
#' 101(4), 882-892. \doi{10.1016/j.jmva.2009.10.009}
#'
#' Besag, J. (1986). On the statistical analysis of dirty pictures. *Journal
#' of the Royal Statistical Society: Series B*, 48(3), 259-279.
#' \doi{10.1111/j.2517-6161.1986.tb01412.x}
#'
#' @export
#' @examples
#' # two regimes on a 10 x 10 lattice: left and right halves
#' set.seed(1)
#' g <- expand.grid(col = 1:10, row = 1:10)
#' D <- nrow(g)
#' W <- 1 * (as.matrix(dist(g)) == 1)
#' dat <- data.frame(x = rnorm(D), vardir = runif(D, 0.3, 0.6))
#' regime <- ifelse(g$col <= 5, 1, 2)
#' dat$y <- ifelse(regime == 1, 2 + 2 * dat$x, 6 - 1 * dat$x) +
#'   rnorm(D, sd = 0.5) + rnorm(D, sd = sqrt(dat$vardir))
#' fit <- scfh_sae(y ~ x, vardir = "vardir", data = dat, W = W, K = 2,
#'                 phi = 0.5, mse = "none", seed = 1)
#' fit
#' table(fit$cluster$cluster, regime)
scfh_sae <- function(
  formula,
  vardir,
  domain = NULL,
  data,
  W,
  K = 2,
  phi = 0.5,
  method = c("AREML", "REML", "ML"),
  n_start = 10,
  maxit = 50,
  mse = c("bootstrap", "none"),
  B = 200,
  seed = NULL,
  conf_level = 0.95
) {
  method <- match.arg(method)
  mse <- match.arg(mse)
  .check_boot_args(mse, B, seed, conf_level)
  .check_count(K, "K")
  .check_count(n_start, "n_start")
  .check_count(maxit, "maxit")
  if (!is.numeric(phi) || length(phi) != 1 || !is.finite(phi) || phi < 0) {
    cli::cli_abort("{.arg phi} must be a single non-negative number.", call = NULL)
  }

  ad <- .prepare_area_data(formula, vardir, domain, data)
  if (length(mgcv::interpret.gam(formula)$smooth.spec) > 0) {
    cli::cli_abort(c(
      "{.arg formula} must not contain smooth terms such as {.code s()}.",
      "i" = "Use {.fn pspline_sae} for spline effects."
    ), call = NULL)
  }
  nb <- .neighbours(W, nrow(data))
  X <- .design_matrix(formula, data)
  p <- ncol(X)
  n_s <- sum(ad$sampled)
  if (K * (p + 1) > n_s) {
    cli::cli_abort(c(
      "{.arg K} = {K} clusters need at least {K * (p + 1)} sampled areas ({p + 1} per cluster), not {n_s}.",
      "i" = "Use fewer clusters or fewer covariates."
    ), call = NULL)
  }
  Xinit <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  Xinit <- Xinit[, apply(Xinit, 2, stats::var) > 0, drop = FALSE]

  fit_once <- function(y) {
    .scfh_fit(y, X, Xinit, ad$vardir, ad$sampled, nb, K, phi, method, n_start, maxit)
  }

  .with_seed(seed, {
    res <- fit_once(ad$y)
    if (!res$admissible) {
      cli::cli_warn(c(
        "The final partition has a cluster with fewer than {p + 1} sampled areas.",
        "i" = "Its parameters are those of an earlier iteration; use a smaller {.arg K}."
      ))
    }
    if (!res$converged) {
      cli::cli_warn("The partition did not stabilize in {maxit} iteration{?s}.")
    }
    boot <- list(mse = rep(NA_real_, nrow(data)), B_success = NA_integer_, n_failed = NA_integer_,
                 errors = NULL, theta = NULL, y = NULL, label = NULL)
    if (mse == "bootstrap") {
      refit <- function(y_star) {
        r <- fit_once(y_star)
        structure(r$est, label = .align_labels(r$cluster, res$cluster, K))
      }
      s2_area <- res$sigma2_v[res$cluster]
      boot <- .pboot_mse(res$synthetic, s2_area, ad$vardir, ad$sampled, refit, B)
    }
  })

  structure(
    list(
      df_est = .est_table(
        domain = ad$domain, direct = ad$y, vardir = ad$vardir, est = res$est,
        synthetic = res$synthetic, random_effect = res$random_effect,
        gamma = res$gamma, mse = boot$mse, sampled = ad$sampled, conf_level = conf_level
      ),
      cluster = data.frame(domain = ad$domain, cluster = res$cluster, sampled = ad$sampled),
      coefficients = res$coefficients,
      sigma2_v = res$sigma2_v,
      size = res$size,
      K = as.integer(K),
      phi = phi,
      loglik = res$loglik,
      objective = res$objective,
      fit = res$fits,
      formula = formula,
      call = match.call(),
      data = data,
      W = W,
      model = "scfh_sae",
      method = method,
      n_start = as.integer(n_start),
      maxit = as.integer(maxit),
      iterations = res$iterations,
      admissible = res$admissible,
      mse_method = mse,
      B = if (mse == "bootstrap") as.integer(B) else NA_integer_,
      B_success = boot$B_success,
      n_boot_failed = boot$n_failed,
      boot_errors = boot$errors,
      boot_theta = boot$theta,
      boot_y = boot$y,
      boot_cluster = boot$label,
      seed = seed,
      conf_level = conf_level,
      convergence = res$converged
    ),
    class = c("flexsae_scfh", "flexsae")
  )
}

#' Neighbour lists from a spatial weight or adjacency matrix
#' @noRd
.neighbours <- function(W, D) {
  if (missing(W) || is.null(W)) {
    cli::cli_abort("{.arg W} is required: a {D} x {D} spatial weight or adjacency matrix.", call = NULL)
  }
  if (inherits(W, "Matrix")) W <- as.matrix(W)
  if (!is.matrix(W) || !is.numeric(W) || nrow(W) != D || ncol(W) != D) {
    cli::cli_abort("{.arg W} must be a numeric {D} x {D} matrix (one row and column per row of {.arg data}).",
                   call = NULL)
  }
  if (anyNA(W) || any(W < 0)) {
    cli::cli_abort("{.arg W} must not contain missing or negative values.", call = NULL)
  }
  A <- (W > 0 | t(W > 0))
  diag(A) <- FALSE
  lapply(seq_len(D), function(d) which(A[d, ]))
}

#' One fit of the SC-FH model on a response vector
#'
#' Used for the point fit and for every bootstrap replicate.
#'
#' @param y Response of all areas (NA for unsampled).
#' @param X Design matrix, all areas.
#' @param Xinit Covariates for the initial k-means (no intercept), all areas.
#' @param nb Neighbour lists over all areas.
#' @return List with `est`, `synthetic`, `random_effect`, `gamma`, `cluster`
#'   (all areas), `sigma2_v`, `size`, `coefficients`, `fits`, `loglik`,
#'   `objective`, `iterations`, `admissible`, `converged`.
#' @noRd
.scfh_fit <- function(y, X, Xinit, vardir, sampled, nb, K, phi, method, n_start, maxit) {
  s <- which(sampled)
  ys <- y[s]
  Xs <- X[s, , drop = FALSE]
  vs <- vardir[s]
  min_size <- ncol(X) + 1 # P + 2 with P covariates (Maranzano et al., 2026, sec. 4.1)
  # neighbours among the sampled areas, in positions of s
  pos <- match(seq_along(sampled), s)
  nb_s <- lapply(nb[s], function(j) stats::na.omit(pos[j]))

  # cluster-wise Fay-Herriot fits; a cluster that is too small or fails keeps
  # its previous parameters (NULL when it never had any)
  fit_clusters <- function(lab, old = vector("list", K)) {
    lapply(seq_len(K), function(k) {
      idx <- which(lab == k)
      if (length(idx) < min_size) {
        return(old[[k]])
      }
      tryCatch(.fit_fh_aliased(ys[idx], Xs[idx, , drop = FALSE], vs[idx], method), error = function(e) old[[k]])
    })
  }
  loglik_matrix <- function(fits) {
    vapply(fits, function(f) {
      if (is.null(f)) {
        return(rep(-Inf, length(s)))
      }
      stats::dnorm(ys, drop(Xs %*% f$beta), sqrt(f$sigma2_v + vs), log = TRUE)
    }, numeric(length(s)))
  }
  objective <- function(lab, ll) {
    same <- sum(vapply(seq_along(lab), function(i) sum(lab[nb_s[[i]]] == lab[i]), numeric(1))) / 2
    sum(ll[cbind(seq_along(lab), lab)]) + phi * same
  }

  if (all(vapply(fit_clusters(rep(1L, length(s))), is.null, logical(1)))) {
    cli::cli_abort("The Fay-Herriot model cannot be fitted to the sampled areas.", call = NULL)
  }
  init <- .scfh_init(if (ncol(Xinit) > 0) Xinit[s, , drop = FALSE] else matrix(ys), K, min_size, n_start)
  lab <- init
  fits <- fit_clusters(lab)
  visited <- list()
  converged <- K == 1
  it <- 0L
  if (K > 1) {
    for (it in seq_len(maxit)) {
      ll <- loglik_matrix(fits)
      new <- lab
      for (i in seq_along(new)) {
        score <- ll[i, ] + phi * tabulate(new[nb_s[[i]]], nbins = K)
        best <- which(score == max(score))
        if (!new[i] %in% best) new[i] <- best[1]
      }
      visited[[length(visited) + 1]] <- list(lab = new, q = objective(new, ll))
      if (identical(new, lab)) {
        converged <- TRUE
        break
      }
      key <- vapply(visited, function(v) paste(v$lab, collapse = ","), character(1))
      if (anyDuplicated(key)) {
        lab <- visited[[which.max(vapply(visited, `[[`, numeric(1), "q"))]]$lab
        converged <- TRUE
        break
      }
      lab <- new
      fits <- fit_clusters(lab, fits)
    }
  }

  # final refit on the final partition; clusters numbered by decreasing size
  size <- tabulate(lab, nbins = K)
  ord <- order(-size, seq_len(K))
  lab <- match(lab, ord)
  fits <- fit_clusters(lab, fits[ord])
  size <- tabulate(lab, nbins = K)
  admissible <- all(size >= min_size) && !any(vapply(fits, is.null, logical(1)))
  ll <- loglik_matrix(fits)
  loglik <- sum(ll[cbind(seq_along(lab), lab)])

  # labels of the unsampled areas
  cluster <- rep(NA_integer_, length(sampled))
  cluster[s] <- lab
  cluster <- .label_unsampled(cluster, nb, Xinit, K)

  fixed <- numeric(length(sampled))
  s2 <- numeric(K)
  for (k in seq_len(K)) {
    idx <- which(cluster == k)
    f <- fits[[k]]
    if (is.null(f)) {
      fixed[idx] <- NA_real_
      next
    }
    fixed[idx] <- drop(X[idx, , drop = FALSE] %*% f$beta)
    s2[k] <- f$sigma2_v
  }
  comb <- .fhrf_combine(fixed, y, s2[cluster], vardir, sampled)

  coefficients <- do.call(rbind, lapply(seq_len(K), function(k) {
    cf <- fits[[k]]$coefficients
    if (is.null(cf)) return(NULL)
    rownames(cf) <- paste0("[", k, "] ", rownames(cf))
    cf
  }))
  list(
    est = comb$est, synthetic = comb$synthetic, random_effect = comb$random_effect,
    gamma = comb$gamma, cluster = cluster, sigma2_v = s2, size = size,
    coefficients = coefficients, fits = fits, loglik = loglik,
    objective = objective(lab, ll),
    iterations = it, admissible = admissible,
    converged = converged && all(vapply(fits, function(f) isTRUE(f$converged), logical(1)))
  )
}

#' Fay-Herriot fit of one cluster without its aliased columns
#'
#' Within a cluster a covariate can be constant or collinear (e.g. a factor
#' level absent from the cluster). As in lm(), such columns are dropped (QR
#' pivoting) and get an NA coefficient and a zero contribution (`beta`).
#' @noRd
.fit_fh_aliased <- function(y, X, vardir, method) {
  q <- qr(X)
  keep <- sort(q$pivot[seq_len(q$rank)])
  f <- .fit_fh_lm(y, X[, keep, drop = FALSE], vardir, method)
  if (length(keep) < ncol(X)) {
    beta <- stats::setNames(numeric(ncol(X)), colnames(X))
    beta[keep] <- f$beta
    cf <- f$coefficients[rep(NA_integer_, ncol(X)), , drop = FALSE]
    cf[keep, ] <- f$coefficients
    rownames(cf) <- colnames(X)
    f$beta <- beta
    f$coefficients <- cf
  }
  f
}

#' Initial partition: k-means on standardized principal components
#'
#' Among `n_start` single-start k-means runs, the admissible one (every
#' cluster with at least `min_size` areas) with the smallest within sum of
#' squares; the best inadmissible one when none is admissible.
#' @noRd
.scfh_init <- function(Z, K, min_size, n_start) {
  n <- nrow(Z)
  if (K == 1) {
    return(rep(1L, n))
  }
  Z <- scale(Z)
  Z[, is.na(colSums(Z))] <- 0
  pc <- stats::prcomp(Z)
  keep <- which(cumsum(pc$sdev^2) / sum(pc$sdev^2) >= 0.9)[1]
  Z <- pc$x[, seq_len(max(keep, 1)), drop = FALSE]
  if (nrow(unique(Z)) < K) {
    cli::cli_abort("The covariates take fewer than {K} distinct values; use a smaller {.arg K}.", call = NULL)
  }
  best <- NULL
  for (i in seq_len(n_start)) {
    km <- stats::kmeans(Z, centers = K, iter.max = 50)
    ok <- all(tabulate(km$cluster, nbins = K) >= min_size)
    if (is.null(best) || (ok && !best$ok) || (ok == best$ok && km$tot.withinss < best$wss)) {
      best <- list(cluster = km$cluster, ok = ok, wss = km$tot.withinss)
    }
  }
  as.integer(best$cluster)
}

#' Labels of the unsampled areas
#'
#' Most frequent label among the labelled neighbours (ties: smallest label),
#' propagated until no label changes; then the nearest cluster mean of the
#' standardized covariates, or cluster 1 (the largest) without covariates.
#' @noRd
.label_unsampled <- function(cluster, nb, Xinit, K) {
  todo <- which(is.na(cluster))
  if (length(todo) == 0) {
    return(cluster)
  }
  repeat {
    changed <- FALSE
    for (d in todo[is.na(cluster[todo])]) {
      lab <- cluster[nb[[d]]]
      lab <- lab[!is.na(lab)]
      if (length(lab) > 0) {
        cluster[d] <- which.max(tabulate(lab, nbins = K))
        changed <- TRUE
      }
    }
    if (!changed || !anyNA(cluster)) break
  }
  left <- which(is.na(cluster))
  if (length(left) > 0) {
    if (ncol(Xinit) == 0) {
      cluster[left] <- 1L
    } else {
      Z <- scale(Xinit)
      Z[, is.na(colSums(Z))] <- 0
      done <- !is.na(cluster)
      centers <- rowsum(Z[done, , drop = FALSE], cluster[done]) / as.vector(table(cluster[done]))
      dist2 <- vapply(seq_len(nrow(centers)), function(k) {
        colSums((t(Z[left, , drop = FALSE]) - centers[k, ])^2)
      }, numeric(length(left)))
      dist2 <- matrix(dist2, nrow = length(left))
      cluster[left] <- as.integer(rownames(centers))[max.col(-dist2, ties.method = "first")]
    }
  }
  as.integer(cluster)
}

#' Map bootstrap labels onto the clusters of the fit (majority rule)
#' @noRd
.align_labels <- function(lab, ref, K) {
  map <- vapply(seq_len(K), function(j) {
    idx <- which(lab == j)
    if (length(idx) == 0) NA_integer_ else which.max(tabulate(ref[idx], nbins = K))
  }, integer(1))
  map[lab]
}

#' Adjusted Rand index (Hubert and Arabie, 1985)
#' @noRd
.adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  c2 <- function(x) sum(x * (x - 1) / 2)
  n <- sum(tab)
  index <- c2(tab)
  expected <- c2(rowSums(tab)) * c2(colSums(tab)) / (n * (n - 1) / 2)
  maximum <- (c2(rowSums(tab)) + c2(colSums(tab))) / 2
  if (maximum == expected) {
    return(1)
  }
  (index - expected) / (maximum - expected)
}
