
# gmethods

<!-- badges: start -->
<!-- badges: end -->

The **gmethods** package implements Robin's g-methods for the estimation of longitudinal causal effects under time-varying confounding. 

## Installation

You can install the development version of `gmethods` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("boyercb/gmethods")
```

## Example

``` r
library(gmethods)

# fit models for the parametric g-formula
gf <- gformula(
  outcome_model = list(
    formula = Y ~ A + L1 + L2 + L3 + lag1_A + lag1_L1 + lag1_L2 + t0,
    link = "logit",
    family = "binomial"
  ),
  covariate_model = list(
    "L1" = list(
      formula = L1 ~ lag1_A + lag1_L1 + lag1_L2 + L3 + t0,
      link = "logit",
      family = "binomial"
    ),
    "L2" = list(
      formula = L2 ~ lag1_A + L1 + lag1_L1 + lag1_L2 + L3 + t0,
      link = "identity",
      family = "normal"
    ),
    "A" = list(
      formula = A ~ lag1_A + L1 + L2 + lag1_L1 + lag1_L2 + L3 + t0,
      link = "logit",
      family = "binomial"
    )
  ),
  data = basicdata_nocomp,
  survival = TRUE,
  id = 'id',
  time = 't0'
)

# simulate population-level interventions
s <- simulate(
  gf,
  interventions = list(
    "Never treat" = function(data, time) set(data, j = "A", value = 0),
    "Always treat" = function(data, time) set(data, j = "A", value = 1)
  ),
  n_samples = 40000
)

# make 'individual' predictions 
 p <- predict(
      gf,
      interventions = list(
        "Never treat" = function(data, time) set(data, j = "A", value = 0),
        "Always treat" = function(data, time) set(data, j = "A", value = 1)
      ),
      n_sims = 500 # number of monte carlo samples per individual
    )

```

