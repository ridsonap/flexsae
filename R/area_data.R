#' Validate and prepare area-level data
#'
#' Shared by the area-level models: parses the formula (mgcv syntax), takes the
#' response, sampling variances and domain labels, and builds the data of the
#' sampled areas used for fitting. Each row of `data` is one area; rows with a
#' missing response are unsampled areas, which are predicted synthetically.
#'
#' @return List with `y`, `vardir`, `domain`, `sampled`, `data` (all areas),
#'   `fit_data` (sampled areas plus internal columns) and `formula_fit`.
#' @noRd
.prepare_area_data <- function(formula, vardir, domain, data) {
  if (!is.data.frame(data) || nrow(data) == 0) {
    cli::cli_abort("{.arg data} must be a data frame with at least one row.", call = NULL)
  }
  if (!inherits(formula, "formula") || length(formula) != 3) {
    cli::cli_abort("{.arg formula} must be a two-sided formula, e.g. {.code y ~ s(x1) + x2}.", call = NULL)
  }
  reserved <- c(.area_var, .resp_var, ".flexsae_w")
  clash <- intersect(reserved, c(names(data), all.vars(formula)))
  if (length(clash) > 0) {
    cli::cli_abort("{.val {clash}} {?is/are} reserved for internal use; rename {?that column/those columns}.", call = NULL)
  }
  gp <- tryCatch(
    mgcv::interpret.gam(formula),
    error = function(e) cli::cli_abort(c("{.arg formula} could not be parsed.", "x" = conditionMessage(e)), call = NULL)
  )
  vars <- all.vars(gp$fake.formula)
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    cli::cli_abort("Variable{?s} {.val {missing_vars}} in {.arg formula} {?is/are} not found in {.arg data}.", call = NULL)
  }

  y <- eval(formula[[2]], data, environment(formula))
  if (!is.numeric(y) || length(y) != nrow(data)) {
    cli::cli_abort("The response must be a numeric vector with one value per row of {.arg data}.", call = NULL)
  }
  y <- as.numeric(y)
  sampled <- !is.na(y)
  if (sum(sampled) < 3) {
    cli::cli_abort("At least 3 sampled areas (non-missing response) are needed, not {sum(sampled)}.", call = NULL)
  }
  if (any(!is.finite(y[sampled]))) {
    cli::cli_abort("The response contains infinite values.", call = NULL)
  }

  rhs_vars <- all.vars(gp$fake.formula[-2])
  na_rows <- which(!stats::complete.cases(data[, rhs_vars, drop = FALSE]))
  if (length(na_rows) > 0) {
    cli::cli_abort(c(
      "Covariates must not contain missing values: {length(na_rows)} row{?s} {?has/have} NA.",
      "i" = "First row{?s}: {.val {utils::head(na_rows, 5)}}."
    ), call = NULL)
  }

  vardir <- .get_variable(data, vardir, "vardir")
  if (!is.numeric(vardir)) {
    cli::cli_abort("{.arg vardir} must be numeric.", call = NULL)
  }
  bad <- which(sampled & (is.na(vardir) | !is.finite(vardir) | vardir <= 0))
  if (length(bad) > 0) {
    cli::cli_abort(c(
      "{.arg vardir} must be positive and finite for every sampled area: {length(bad)} area{?s} {?is/are} not.",
      "i" = "First row{?s}: {.val {utils::head(bad, 5)}}."
    ), call = NULL)
  }

  if (is.null(domain)) {
    domain <- seq_len(nrow(data))
  } else {
    domain <- .get_variable(data, domain, "domain")
    dup <- unique(domain[duplicated(domain)])
    if (length(dup) > 0) {
      cli::cli_warn(c(
        "{.arg domain} has {length(dup)} duplicated value{?s}: {.val {utils::head(dup, 10)}}.",
        "i" = "Each row is fitted as a separate area; the domain column only labels the results."
      ))
    }
  }

  fit_data <- data[sampled, , drop = FALSE]
  fit_data[[.resp_var]] <- y[sampled]
  fit_data[[.area_var]] <- factor(which(sampled))
  fit_data$.flexsae_w <- 1 / vardir[sampled]

  list(
    y = y,
    vardir = as.numeric(vardir),
    domain = domain,
    sampled = sampled,
    data = data,
    fit_data = fit_data,
    formula_fit = .fh_formula(formula)
  )
}
