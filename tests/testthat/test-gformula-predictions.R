
# survival ----------------------------------------------------------------

gf <- gformula(
  outcome_model = list(
    formula = Y ~ A + L1 + L2 + L3 + lag1_A + lag1_L1 + lag1_L2 + t0,
    link = "logit",
    family = "binomial"
  ),
  covariate_model = list(
    "L1" = list(
      formula = L1 ~ lag1_A + lag1_cumavg_L1 + lag1_cumavg_L2 + L3 + t0,
      link = "logit",
      family = "binomial"
    ),
    "L2" = list(
      formula = L2 ~ lag1_A + L1 + lag1_cumavg_L1 + lag1_cumavg_L2 + L3 + t0,
      link = "identity",
      family = "normal"
    ),
    "A" = list(
      formula = A ~ lag1_A + L1 + L2 + lag1_cumavg_L1 + lag1_cumavg_L2 + L3 + t0,
      link = "logit",
      family = "binomial"
    )
  ),
  data = basicdata_nocomp,
  survival = TRUE,
  id = 'id',
  time = 't0'
)

test_that("Prediction arguments work as expected.", {
  # basic usage (default sims, data same as fit)
  expect_silent({
    p <- predict(gf)

    p <- predict(
      object = gf,
      interventions = list(
        "Never treat" = function(data, time) set(data, j = "A", value = 0),
        "Always treat" = function(data, time) set(data, j = "A", value = 1)
      )
    )

    p <- predict(
      object = gf,
      interventions = list(
        "Never treat" = function(data, time) set(data, j = "A", value = 0),
        "Always treat" = function(data, time) set(data, j = "A", value = 1)
      ),
      last_only = FALSE
    )

    p <- predict(
      object = gf,
      interventions = list(
        "Never treat" = function(data, time) set(data, j = "A", value = 0),
        "Always treat" = function(data, time) set(data, j = "A", value = 1)
      ),
      last_only = FALSE,
      predict_covs = TRUE
    )

    p <- predict(
      object = gf,
      interventions = list(
        "Never treat" = function(data, time) set(data, j = "A", value = 0),
        "Always treat" = function(data, time) set(data, j = "A", value = 1)
      ),
      n_boots = 10,
      last_only = FALSE,
      predict_covs = TRUE
    )
  })
})




# binary eof --------------------------------------------------------------



# continuous eof ----------------------------------------------------------


