#' Compute B-spline coefficients using ridge regression
#'
#' @param data A data.frame with columns: ID, Time, Value
#' @param knots Numeric vector of internal knot positions for splines::bs()
#' @param degree Integer spline degree. Default is 3 (cubic).
#' @param lambda Numeric lambda for glmnet (ridge). Default is 1e-7.
#' @param standardize Logical; passed to glmnet::glmnet(). Default TRUE.
#'
#' @return A named numeric vector of coefficients: (Intercept) and bs_1..bs_p.
#'         If fitting is not possible, returns an NA vector of the same length.
#' @importFrom splines bs
#' @importFrom glmnet glmnet
#' @importFrom stats setNames coef
#' @keywords internal
calculate_bspline_coefficients <- function(data, knots, degree = 3, lambda = 1e-7, standardize = TRUE) {
  stopifnot(all(c("Time","Value") %in% names(data)))

  x <- tryCatch(
    splines::bs(data$Time, knots = knots, degree = degree, intercept = FALSE),
    error = function(e) NULL
  )

  if (is.null(x)) {
    nm <- c("(Intercept)")
    return(stats::setNames(rep(NA_real_, length(nm)), nm))
  }

  x <- as.matrix(x)
  p <- ncol(x)
  nm <- c("(Intercept)", if (p > 0) paste0("bs_", seq_len(p)) else character())

  if (p == 0 || nrow(x) < 1L || length(data$Value) != nrow(x)) {
    return(stats::setNames(rep(NA_real_, length(nm)), nm))
  }

  y <- data$Value

  fit <- tryCatch(
    glmnet::glmnet(x, y, alpha = 0, lambda = lambda, standardize = standardize, intercept = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(stats::setNames(rep(NA_real_, length(nm)), nm))
  }

  cf <- as.vector(stats::coef(fit))

  if (length(cf) < length(nm)) {
    cf <- c(cf, rep(NA_real_, length(nm) - length(cf)))
  } else if (length(cf) > length(nm)) {
    cf <- cf[seq_along(nm)]
  }
  names(cf) <- nm
  cf
}
