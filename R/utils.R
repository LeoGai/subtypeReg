#' Safe lm.fit wrapper
#' @keywords internal
.lm_fit_safe <- function(X, y) {
  z <- stats::lm.fit(X, y)
  z$coefficients[is.na(z$coefficients)] <- 0
  z
}

#' Expand and sanitize candidate k values
#' @keywords internal
.valid_k_range <- function(k_range, n_id) {
  k_range <- sort(unique(k_range))
  k_range <- k_range[k_range >= 2 & k_range <= max(2, n_id)]
  if (length(k_range) == 0) k_range <- 2:min(8, n_id)
  k_range
}
