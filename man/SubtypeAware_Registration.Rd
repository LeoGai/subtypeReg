\name{SubtypeAware_Registration}
\alias{SubtypeAware_Registration}
\title{Subtype-aware Registration of Longitudinal Data}
\description{
Registers subject-specific time axes by fitting spline coefficients under multiple temporal shifts,
clustering baseline coefficients, selecting representative subjects, and choosing the best shift per subject
to align trajectories.
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
  lambda = 1e-07,
  verbose = FALSE
)
}
\arguments{
  \item{data}{A \code{data.frame} with columns \code{ID}, \code{Feature}, \code{Time}, \code{Value}.}
  \item{alpha}{Numeric in (0, 1]; the quantile used to select representative points within each cluster.}
  \item{k_range}{Integer vector of candidate cluster numbers, e.g., \code{2:8}.}
  \item{timepos_option}{Numeric vector of temporal shifts to consider (can be negative), e.g., \code{c(-1, 0, 1)} or \code{c(4, 3, 2, 1)}. The order determines the tie-break if distances are equal.}
  \item{tau}{Numeric threshold on the initial average silhouette width; if the initial value exceeds \code{tau}, fewer iterations are performed.}
  \item{tmin}{Numeric, lower bound of the time window.}
  \item{tmax}{Numeric, upper bound of the time window.}
  \item{knots}{Numeric vector of internal knot positions passed to \code{splines::bs}.}
  \item{degree}{Integer spline degree passed to \code{splines::bs}; default is \code{3}.}
  \item{lambda}{Numeric ridge penalty passed to \code{glmnet}; default is \code{1e-7}.}
  \item{verbose}{Logical; if \code{TRUE}, prints silhouette vectors at each iteration and the final
  selected \code{k} and cluster centers (means).}
}
\value{
A \code{data.frame} of the same structure as \code{data}, with \code{Time} adjusted (registered) and
restricted to the interval \code{[tmin, tmax]}.
}
\details{
For each subject (\code{ID}), spline coefficients are estimated at \code{Original} and at each
candidate shift in \code{timepos_option}. Baseline coefficients (\code{Original}) are clustered via
\code{cluster::pam}; representative subjects per cluster are selected by a quantile cutoff \code{alpha}.
For each selected subject, the best state (shift) is chosen by minimizing squared distance to the local
cluster center; ties are broken by the order of \code{timepos_option} (with \code{Original} first).
An optional iterative refinement updates clusters and selections until the stopping rule is met:
no change in \code{k} and no improvement in the top-two silhouette scores, or the iteration limit is reached.
}
\examples{
set.seed(1)
toy <- do.call(rbind, lapply(1:6, function(id) {
  t <- 0:10
  data.frame(
    ID = id, Feature = "X", Time = t,
    Value = sin(t/3) + rnorm(length(t), sd = 0.1)
  )
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
  lambda = 1e-7,
  verbose = FALSE
)

head(out)
}

