#' Compute B-spline coefficients (ridge by glmnet, with fallback)
#'
#' @param data A data.frame with columns: ID, Time, Value
#' @param knots Numeric vector of internal knot positions for splines::bs
#' @param degree Integer spline degree. Default 3 (cubic).
#' @param lambda Numeric lambda for glmnet (ridge). Default 1e-7.
#'
#' @return Named numeric vector of coefficients: (Intercept) + bs_* columns.
#' @importFrom splines bs
#' @importFrom glmnet glmnet
#' @keywords internal
calculate_bspline_coefficients <- function(data, knots, degree = 3, lambda = 1e-7) {
  stopifnot(all(c("Time","Value") %in% names(data)))
  x <- tryCatch(
    splines::bs(data$Time, knots = knots, degree = degree, intercept = FALSE),
    error = function(e) NULL
  )
  if (is.null(x) || nrow(as.matrix(x)) == 0 ||
      length(data$Value) < (ncol(as.matrix(x)) + 1L)) {
    p <- if (!is.null(x)) ncol(as.matrix(x)) else 0L
    nm <- c("(Intercept)", if (p > 0) paste0("bs_", seq_len(p)) else character())
    out <- rep(NA_real_, length(nm)); names(out) <- nm
    return(out)
  }
  x <- as.matrix(x)
  y <- data$Value

  fit <- tryCatch(
    glmnet::glmnet(x, y, alpha = 0, lambda = lambda, standardize = FALSE),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    cf <- as.vector(glmnet::coef.glmnet(fit))
    p <- ncol(x)
    nm <- c("(Intercept)", paste0("bs_", seq_len(p)))
    names(cf) <- nm
    return(cf)
  } else {
    X <- cbind(Intercept = 1, x)
    if (nrow(X) <= ncol(X)) {
      p <- ncol(X)
      nm <- colnames(X)
      out <- rep(NA_real_, p); names(out) <- nm
      return(out)
    }
    cf <- tryCatch(
      {
        fit_lm <- .lm_fit_safe(X, y)
        as.numeric(fit_lm$coefficients)
      },
      error = function(e) rep(NA_real_, ncol(X))
    )
    names(cf) <- colnames(X)
    return(cf)
  }
}
