#' @importFrom stats rnorm
#' @importFrom utils head
NULL

utils::globalVariables(c(".flexsae_w"))

`%||%` <- function(x, y) if (is.null(x)) y else x

# Name of the internal area factor and response column added to the fitting data
.area_var <- ".flexsae_area"
.resp_var <- ".flexsae_y"

#' Get a variable from data by column name, one-sided formula or vector
#'
#' Same contract as `.get_variable()` in fastsae.
#' @noRd
.get_variable <- function(data, variable, arg = "variable") {
  if (is.character(variable) && length(variable) == 1) {
    if (!variable %in% colnames(data)) {
      cli::cli_abort("Column {.val {variable}} given in {.arg {arg}} is not found in {.arg data}.")
    }
    return(data[[variable]])
  }
  if (inherits(variable, "formula")) {
    v_names <- all.vars(variable)
    if (length(v_names) != 1 || !v_names %in% colnames(data)) {
      cli::cli_abort("{.arg {arg}} must be a one-sided formula naming one column of {.arg data}.")
    }
    return(data[[v_names]])
  }
  if (length(variable) == nrow(data)) {
    return(variable)
  }
  cli::cli_abort(
    "{.arg {arg}} must be a column name, a one-sided formula or a vector of length {nrow(data)}, not length {length(variable)}."
  )
}

#' Evaluate code with a fixed seed, then restore the random number stream
#'
#' A NULL seed leaves the stream as it is.
#' @noRd
.with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(code)
  }
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  })
  set.seed(seed)
  code
}

#' Add the area random effect term to the right-hand side of a formula
#'
#' The response is replaced by the internal response column, so that bootstrap
#' replicates only swap one column of the fitting data.
#' @noRd
.fh_formula <- function(formula) {
  f <- formula
  f[[2]] <- as.name(.resp_var)
  f[[3]] <- call("+", f[[3]], call("s", as.name(.area_var), bs = "re"))
  f
}

#' Standard table of area estimates
#'
#' Every model fills the same columns; the ones that do not apply are NA.
#' @noRd
.est_table <- function(domain, direct, vardir, est, synthetic, random_effect,
                       gamma, mse, sampled, conf_level) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  data.frame(
    domain = domain,
    direct = direct,
    vardir = vardir,
    est = est,
    synthetic = synthetic,
    random_effect = random_effect,
    gamma = gamma,
    mse = mse,
    rse = 100 * sqrt(mse) / abs(est),
    ci_lower = est - z * sqrt(mse),
    ci_upper = est + z * sqrt(mse),
    sampled = sampled,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
