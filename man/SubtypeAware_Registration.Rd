\
\name{SubtypeAware_Registration}
\alias{SubtypeAware_Registration}
\title{Subtype-aware temporal registration via spline-coefficient clustering}
\description{
Perform temporal registration by fitting B-splines under multiple time shifts,
embedding each subject by spline coefficients, clustering with PAM, and selecting
the best shift per subject.
}
\usage{
SubtypeAware_Registration(
  data,
  alpha = 0.95,
  k_range = 2:8,
  timepos_option = c(4, 3, 2, 1),
  tau = 0.45,
  tmin,
  tmax,
  knots,
  degree = 3,
  lambda = 1e-7
)
}
\arguments{
\item{data}{A \code{data.frame} with columns \code{ID}, \code{Feature}, \code{Time}, \code{Value}.}

\item{alpha}{Quantile in (0,1] to select representative points per cluster.}

\item{k_range}{Integer vector of candidate cluster numbers, e.g., \code{2:8}.}

\item{timepos_option}{Numeric vector of candidate time shifts (can be negative).}

\item{tau}{Silhouette threshold. If the initial mean silhouette exceeds \code{tau}, fewer iterations are performed.}

\item{tmin}{Lower bound of the time window.}

\item{tmax}{Upper bound of the time window.}

\item{knots}{Numeric vector of internal knot positions for \code{splines::bs}.}

\item{degree}{Spline degree (default 3).}

\item{lambda}{Ridge penalty for \code{glmnet} (default \code{1e-7}).}
}
\value{
A \code{data.frame} with the same columns as \code{data}, where \code{Time}
has been shifted by the selected optimal amount per subject and truncated to \code{[tmin, tmax]}.
}
\details{
Only the exported user-level function is documented here. Helper functions are internal.
}
\examples{
set.seed(1)
toy <- do.call(rbind, lapply(1:6, function(id) {
  t <- 0:10
  data.frame(ID = id, Feature = "X", Time = t,
             Value = sin(t/3) + rnorm(length(t), sd = 0.1))
}))

out <- SubtypeAware_Registration(
  data = toy,
  alpha = 0.95,
  k_range = 2:4,
  timepos_option = c(-1, 0, 1),
  tau = 0.45,
  tmin = 0, tmax = 10,
  knots = c(3, 7),
  degree = 3,
  lambda = 1e-7
)
head(out)

\donttest{
out2 <- SubtypeAware_Registration(
  data = toy, alpha = 0.9, k_range = 2:5, timepos_option = c(-2,-1,0,1,2),
  tau = 0.45, tmin = 0, tmax = 10, knots = c(3,7), degree = 3, lambda = 1e-7
)
}
}
