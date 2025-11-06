# subtypeReg

Subtype-Aware Registration of Longitudinal Electronic Health Records.

## Installation

```r
# install.packages("devtools")
devtools::install_github("LeoGai/subtypeReg")
```

## Example

```r
library(subtypeReg)

set.seed(1)
toy <- do.call(rbind, lapply(1:6, function(id) {
  t <- 0:10
  data.frame(ID = id, Feature = "X", Time = t, Value = sin(t/3) + rnorm(length(t), sd = 0.1))
}))

out <- SubtypeAware_Registration(
  data = toy,
  alpha = 0.95,
  k_range = 2:4,
  timepos_option = c(-1,0,1),
  tau = 0.45,
  tmin = 0, tmax = 10,
  knots = c(3,7),
  degree = 3,
  lambda = 1e-7
)

head(out)
```
