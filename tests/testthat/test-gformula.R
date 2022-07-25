# if (!requireNamespace("gfoRmula", quietly = TRUE, exclude)) {
#   stop("This test uses \"gfoRmula\" package.")
# }

library(gfoRmula, exclude = c("gformula"), warn.conflicts = FALSE)

# scenario 1 --------------------------------------------------------------

## Estimating the effect of static treatment strategies on risk of a
## failure event

# # gfoRmula --------------------------------------------------------------

id <- 'id'
time_points <- 7
time_name <- 't0'
covnames <- c('L1', 'L2', 'A')
outcome_name <- 'Y'
covtypes <- c('binary', 'normal', 'binary')
histories <- c(lagged)
histvars <- list(c('A', 'L1', 'L2'))
covparams <- list(covmodels = c(L1 ~ lag1_A + lag1_L1 + lag1_L2 +
                                  L3 + t0,
                                L2 ~ lag1_A + L1 + lag1_L1 + lag1_L2
                                + L3 + t0,
                                A ~ lag1_A + L1 + L2 + lag1_L1 +
                                  lag1_L2 + L3 + t0))
ymodel <- Y ~ A + L1 + L2 + L3 + lag1_A + lag1_L1 + lag1_L2 + t0
intvars <- list('A', 'A')
interventions <- list(list(c(static, rep(0, time_points))),
                      list(c(static, rep(1, time_points))))
int_descript <- c('Never treat', 'Always treat')
nsimul <- 40000

gform_basic_nocomp <- gformula_survival(obs_data = basicdata_nocomp, id = id,
                                        time_points = time_points,
                                        time_name = time_name, covnames = covnames,
                                        outcome_name = outcome_name,
                                        covtypes = covtypes,
                                        covparams = covparams, ymodel = ymodel,
                                        intvars = intvars,
                                        interventions = interventions,
                                        int_descript = int_descript,
                                        histories = histories, histvars = histvars,
                                        basecovs = c('L3'), nsimul = nsimul,
                                        model_fits = TRUE,
                                        seed = 1234)


# # custom ----------------------------------------------------------------

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

s <- simulate(
  gf,
  interventions = list(
    "Never treat" = function(data, time) set(data, j = "A", value = 0),
    "Always treat" = function(data, time) set(data, j = "A", value = 1)
  ),
  n_samples = 40000
)


# # test ------------------------------------------------------------------

# check that model coefficients are the same
test_that("scenario 1: models fits are correct", {
  # outcome
  expect_equal(gf$fit$outcome$coefficients,
               gform_basic_nocomp$fits$Y$coefficients)

  # covariates
  expect_equal(gf$fit$covariate$L1$coefficients,
               gform_basic_nocomp$fits$L1$coefficients)

  expect_equal(gf$fit$covariate$L2$coefficients,
               gform_basic_nocomp$fits$L2$coefficients)

  expect_equal(gf$fit$covariate$A$coefficients,
               gform_basic_nocomp$fits$A$coefficients)

})


test_that("scenario 1: natural course is the same", {
  nprisk_new <- s$np_means$obs_means
  nprisk_old <- gform_basic_nocomp$result

  expect_equal(
    nprisk_new$poprisk[nprisk_new$t0 == 6],
    nprisk_old$`NP Risk`[nprisk_old$k == 6 & nprisk_old$Interv. == 0]
  )
})

test_that("scenario 1: results are close", {
  result_new <- s$results
  result_old <- gform_basic_nocomp$result[gform_basic_nocomp$result$k == 6, ]

  # risks
  expect_equal(
    result_new$means$estimate,
    result_old$`g-form risk`,
    tolerance = 0.02
  )

  result_new$ratios$estimate[is.na(result_new$ratios$estimate)] <- 1

  # ratios
  expect_equal(
    result_new$ratios$estimate,
    result_old$`Risk ratio`,
    tolerance = 0.02
  )

  result_new$diffs$estimate[is.na(result_new$diffs$estimate)] <- 0

  # differences
  expect_equal(
    result_new$diffs$estimate,
    result_old$`Risk difference`,
    tolerance = 0.02
  )
})


# scenario 2 --------------------------------------------------------------

## Estimating the effect of treatment strategies on risk of a failure event
## when competing events exist

# # gfoRmula --------------------------------------------------------------

id <- 'id'
time_points <- 7
time_name <- 't0'
covnames <- c('L1', 'L2', 'A')
outcome_name <- 'Y'
compevent_name <- 'D'
covtypes <- c('binary', 'normal', 'binary')
histories <- c(lagged)
histvars <- list(c('A', 'L1', 'L2'))
covparams <- list(covlink = c('logit', 'identity', 'logit'),
                  covmodels = c(L1 ~ lag1_A + lag1_L1 + lag1_L2 +
                                  L3 + t0,
                                L2 ~ lag1_A + L1 + lag1_L1 + lag1_L2
                                  + L3 + t0,
                                A ~ lag1_A + L1 + L2 + lag1_L1 +
                                  lag1_L2 + L3 + t0))
ymodel <- Y ~ A + L1 + L2 + lag1_A + lag1_L1 + lag1_L2 + L3 + t0
compevent_model <- D ~ A + L1 + L2 + lag1_A + lag1_L1 + lag1_L2 + L3 + t0
intvars <- list('A', 'A')
interventions <- list(list(c(static, rep(0, time_points))),
                      list(c(static, rep(1, time_points))))
int_descript <- c('Never treat', 'Always treat')
nsimul <- 40000

gform_basic <- gformula_survival(obs_data = basicdata, id = id,
                                 time_points = time_points,
                                 time_name = time_name, covnames = covnames,
                                 outcome_name = outcome_name,
                                 compevent_name = compevent_name,
                                 covtypes = covtypes,
                                 covparams = covparams, ymodel = ymodel,
                                 compevent_model = compevent_model,
                                 intvars = intvars, interventions = interventions,
                                 int_descript = int_descript,
                                 histories = histories, histvars = histvars,
                                 model_fits = TRUE,
                                 basecovs = c('L3'), nsimul = nsimul,
                                 seed = 1234)


# # custom ----------------------------------------------------------------

gf <- gformula(
  outcome_model = list(
    formula = Y ~ A + L1 + L2 + lag1_A + lag1_L1 + lag1_L2 + L3 + t0,
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
  compevent_model = list(
    formula = D ~ A + L1 + L2 + lag1_A + lag1_L1 + lag1_L2 + L3 + t0,
    link = "logit",
    family = "binomial"
  ),
  censor_compevent = FALSE,
  data = basicdata,
  survival = TRUE,
  id = 'id',
  time = 't0'
)

s <- simulate(
  gf,
  interventions = list(
    "Never treat" = function(data, time) set(data, j = "A", value = 0),
    "Always treat" = function(data, time) set(data, j = "A", value = 1)
  ),
  n_samples = 40000
)


# # test ------------------------------------------------------------------

# check that model coefficients are the same
test_that("scenario 2: models fits are correct", {
  # outcome
  expect_equal(gf$fit$outcome$coefficients,
               gform_basic$fits$Y$coefficients)

  # covariates
  expect_equal(gf$fit$covariate$L1$coefficients,
               gform_basic$fits$L1$coefficients)

  expect_equal(gf$fit$covariate$L2$coefficients,
               gform_basic$fits$L2$coefficients)

  expect_equal(gf$fit$covariate$A$coefficients,
               gform_basic$fits$A$coefficients)

})


test_that("scenario 2: natural course is the same", {
  nprisk_new <- s$np_means$obs_means
  nprisk_old <- gform_basic$result

  expect_equal(
    nprisk_new$poprisk[nprisk_new$t0 == 6],
    nprisk_old$`NP Risk`[nprisk_old$k == 6 & nprisk_old$Interv. == 0]
  )
})

test_that("scenario 2: results are close", {
  result_new <- s$results
  result_old <- gform_basic$result[gform_basic$result$k == 6, ]

  # risks
  expect_equal(
    result_new$means$estimate,
    result_old$`g-form risk`,
    tolerance = 0.02
  )

  result_new$ratios$estimate[is.na(result_new$ratios$estimate)] <- 1

  # ratios
  expect_equal(
    result_new$ratios$estimate,
    result_old$`Risk ratio`,
    tolerance = 0.02
  )

  result_new$diffs$estimate[is.na(result_new$diffs$estimate)] <- 0

  # differences
  expect_equal(
    result_new$diffs$estimate,
    result_old$`Risk difference`,
    tolerance = 0.02
  )
})


# scenario 3 --------------------------------------------------------------

## Using IP weighting to estimate natural course risk
## Only the natural course intervention is included for simplicity

# # gfoRmula --------------------------------------------------------------

covnames <- c('L', 'A')
histories <- c(lagged)
histvars <- list(c('A', 'L'))
ymodel <- Y ~ L + A
covtypes <- c('binary', 'normal')
covparams <- list(covmodels = c(L ~ lag1_L + lag1_A,
                                A ~ lag1_L + L + lag1_A))
censor_name <- 'C'
censor_model <- C ~ L
res_censor <- gformula_survival(obs_data = censor_data, id = 'id',
                                time_name = 't0', covnames = covnames,
                                outcome_name = 'Y',
                                censor_name = censor_name, censor_model = censor_model,
                                covtypes = covtypes,
                                covparams = covparams, ymodel = ymodel,
                                intvars = NULL, interventions = NULL,
                                int_descript = NULL, model_fits = TRUE,
                                histories = histories, histvars = histvars,
                                seed = 1234)


# # custom ----------------------------------------------------------------

gf <- gformula(
  outcome_model = list(
    formula = Y ~ L + A,
    link = "logit",
    family = "binomial"
  ),
  covariate_model = list(
    "L" = list(
      formula = L ~ lag1_L + lag1_A,
      link = "logit",
      family = "binomial"
    ),
    "A" = list(
      formula = A ~ lag1_L + L + lag1_A,
      link = "identity",
      family = "normal"
    )
  ),
  censor_model = list(
    formula = C ~ L,
    link = "logit",
    family = "binomial"
  ),
  data = censor_data,
  survival = TRUE,
  id = 'id',
  time = 't0'
)

s <- simulate(
  gf,
  interventions = NULL
)


# # test ------------------------------------------------------------------

# check that model coefficients are the same
test_that("scenario 3: models fits are correct", {
  # outcome
  expect_equal(gf$fit$outcome$coefficients,
               res_censor$fits$Y$coefficients)

  # covariates
  expect_equal(gf$fit$covariate$L1$coefficients,
               res_censor$fits$L1$coefficients)

  expect_equal(gf$fit$covariate$L2$coefficients,
               res_censor$fits$L2$coefficients)

  expect_equal(gf$fit$covariate$A$coefficients,
               res_censor$fits$A$coefficients)

})


test_that("scenario 3: natural course is the same", {
  nprisk_new <- s$np_means$obs_means
  nprisk_old <- res_censor$result

  expect_equal(
    nprisk_new$poprisk[nprisk_new$t0 == 9],
    nprisk_old$`IP weighted risk`[nprisk_old$k == 9 & nprisk_old$Interv. == 0]
  )
})

test_that("scenario 3: results are close", {
  result_new <- s$results
  result_old <- res_censor$result[res_censor$result$k == 9, ]

  # risks
  expect_equal(
    result_new$means$estimate,
    result_old$`IP weighted risk`,
    tolerance = 0.02
  )

  result_new$ratios$estimate[is.na(result_new$ratios$estimate)] <- 1

  # ratios
  expect_equal(
    result_new$ratios$estimate,
    result_old$`Risk ratio`,
    tolerance = 0.02
  )

  result_new$diffs$estimate[is.na(result_new$diffs$estimate)] <- 0

  # differences
  expect_equal(
    result_new$diffs$estimate,
    result_old$`Risk difference`,
    tolerance = 0.02
  )
})


# scenario 4: -------------------------------------------------------------

## Estimating the effect of treatment strategies on the mean of a continuous
## end of follow-up outcome


# # gfoRmula --------------------------------------------------------------

id <- 'id'
time_name <- 't0'
covnames <- c('L1', 'L2', 'A')
outcome_name <- 'Y'
covtypes <- c('categorical', 'normal', 'binary')
histories <- c(lagged)
histvars <- list(c('A', 'L1', 'L2'))
covparams <- list(covmodels = c(L1 ~ lag1_A + lag1_L1 + L3 + t0 +
                                 lag1_L2,
                                L2 ~ lag1_A + L1 + lag1_L1 + lag1_L2 + L3 + t0,
                                A ~ lag1_A + L1 + L2 + lag1_L1 + lag1_L2 + L3 + t0))
ymodel <- Y ~ A + L1 + L2 + lag1_A + lag1_L1 + lag1_L2 + L3
intvars <- list('A', 'A')
interventions <- list(list(c(static, rep(0, 7))),
                      list(c(static, rep(1, 7))))
int_descript <- c('Never treat', 'Always treat')
nsimul <- 50000

gform_cont_eof <- gformula_continuous_eof(obs_data = continuous_eofdata,
                                          id = id,
                                          time_name = time_name,
                                          covnames = covnames,
                                          outcome_name = outcome_name,
                                          covtypes = covtypes,
                                          covparams = covparams, ymodel = ymodel,
                                          intvars = intvars,
                                          interventions = interventions,
                                          int_descript = int_descript,
                                          histories = histories,
                                          histvars = histvars,
                                          model_fits = TRUE,
                                          basecovs = c("L3"),
                                          nsimul = nsimul, seed = 1234)


# # custom ----------------------------------------------------------------

gf <- gformula(
  outcome_model = list(
    formula = Y ~ A + L1 + L2 + lag1_A + lag1_L1 + lag1_L2 + L3,
    link = "identity",
    family = "normal"
  ),
  covariate_model = list(
    "L1" = list(
      formula = L1 ~ lag1_A + lag1_L1 + L3 + t0 + lag1_L2,
      family = "categorical"
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
  data = continuous_eofdata,
  survival = FALSE,
  id = 'id',
  time = 't0'
)

s <- simulate(
  gf,
  interventions = list(
    "Never treat" = function(data, time) set(data, j = "A", value = 0),
    "Always treat" = function(data, time) set(data, j = "A", value = 1)
  ),
  n_samples = 50000
)


# # test ------------------------------------------------------------------

# check that model coefficients are the same
test_that("scenario 4: models fits are correct", {
  # outcome
  expect_equal(gf$fit$outcome$coefficients,
               gform_cont_eof$fits$Y$coefficients)

  # covariates
  expect_equal(gf$fit$covariate$L1$coefficients,
               gform_cont_eof$fits$L1$coefficients)

  expect_equal(gf$fit$covariate$L2$coefficients,
               gform_cont_eof$fits$L2$coefficients)

  expect_equal(gf$fit$covariate$A$coefficients,
               gform_cont_eof$fits$A$coefficients)

})


test_that("scenario 4: natural course is the same", {
  nprisk_new <- s$np_means$obs_means
  nprisk_old <- gform_cont_eof$result

  expect_equal(
    nprisk_new$Ey,
    nprisk_old$`NP mean`[nprisk_old$k == 6 & nprisk_old$Interv. == 0]
  )
})

test_that("scenario 4: results are close", {
  result_new <- s$results
  result_old <- gform_cont_eof$result[gform_cont_eof$result$k == 6, ]

  # risks
  expect_equal(
    result_new$means$estimate,
    result_old$`g-form mean`,
    tolerance = 0.02
  )

  result_new$ratios$estimate[is.na(result_new$ratios$estimate)] <- 1

  # ratios
  expect_equal(
    result_new$ratios$estimate,
    result_old$`Mean ratio`,
    tolerance = 0.02
  )

  result_new$diffs$estimate[is.na(result_new$diffs$estimate)] <- 0

  # differences
  expect_equal(
    result_new$diffs$estimate,
    result_old$`Mean difference`,
    tolerance = 0.03
  )
})


# scenario 5: -------------------------------------------------------------

## Estimating the effect of threshold interventions on the mean of a binary
## end of follow-up outcome

# # gfoRmula --------------------------------------------------------------

id <- 'id_num'
time_name <- 'time'
covnames <- c('cov1', 'cov2', 'treat')
outcome_name <- 'outcome'
histories <- c(lagged, cumavg)
histvars <- list(c('treat', 'cov1', 'cov2'), c('cov1', 'cov2'))
covtypes <- c('binary', 'zero-inflated normal', 'normal')
covparams <- list(covmodels = c(cov1 ~ lag1_treat + lag1_cov1 + lag1_cov2 + cov3 +
                                  time,
                                cov2 ~ lag1_treat + cov1 + lag1_cov1 + lag1_cov2 +
                                  cov3 + time,
                                treat ~ lag1_treat + cumavg_cov1 +
                                  cumavg_cov2 + cov3 + time))
ymodel <- outcome ~  treat + cov1 + cov2 + lag1_cov1 + lag1_cov2 + cov3
intvars <- list('treat', 'treat')
interventions <- list(list(c(static, rep(0, 7))),
                      list(c(threshold, 1, Inf)))
int_descript <- c('Never treat', 'Threshold - lower bound 1')
nsimul <- 10000

gform_bin_eof <- gformula_binary_eof(obs_data = binary_eofdata, id = id,
                                     time_name = time_name,
                                     covnames = covnames,
                                     outcome_name = outcome_name,
                                     covtypes = covtypes,
                                     covparams = covparams,
                                     ymodel = ymodel,
                                     intvars = intvars,
                                     interventions = interventions,
                                     int_descript = int_descript,
                                     histories = histories, histvars = histvars,
                                     basecovs = c("cov3"), seed = 152354,
                                     model_fits = TRUE, sim_data_b = TRUE,
                                     nsimul = 100000)


# # custom ----------------------------------------------------------------

d <- binary_eofdata
d[, t0 := time]
d[, time := NULL]

gf <- gformula(
  outcome_model = list(
    formula = outcome ~ treat + cov1 + cov2 + lag1_cov1 + lag1_cov2 + cov3,
    link = "logit",
    family = "binomial"
  ),
  covariate_model = list(
    "cov1" = list(
      formula = cov1 ~ lag1_treat + lag1_cov1 + lag1_cov2 + cov3 + t0,
      link = "logit",
      family = "binomial"
    ),
    "cov2" = list(
      formula = cov2 ~ lag1_treat + cov1 + lag1_cov1 + lag1_cov2 +
        cov3 + t0,
      link = "identity",
      family = "zero-inflated normal"
    ),
    "treat" = list(
      formula = treat ~ lag1_treat + cumavg_cov1 + cumavg_cov2 +
        cov3 + t0,
      link = "identity",
      family = "normal"
    )
  ),
  data = d,
  survival = FALSE,
  id = 'id_num',
  time = 't0'
)

s <- simulate(
  gf,
  interventions = list(
    "Never treat" = function(data, time) set(data, j = "treat", value = 0),
    "Threshold - lower bound 1" = function(data, time)
      set(data,
          j = "treat",
          value = ifelse(data$treat < 1, 1, data$treat))
  ),
  n_samples = 10000,
  return_sims = TRUE,
  bound_sims = TRUE
)


# # test ------------------------------------------------------------------

# check that model coefficients are the same
test_that("scenario 5: models fits are correct", {

  names(gf$fit$covariate$cov1$coefficients) <-
    gsub("t0", "time", names(gf$fit$covariate$cov1$coefficients))

  names(gf$fit$covariate$cov2$fits$zero$coefficients) <-
    gsub("t0", "time", names(gf$fit$covariate$cov2$fits$zero$coefficients))

  names(gf$fit$covariate$cov2$fits$normal$coefficients) <-
    gsub("t0", "time", names(gf$fit$covariate$cov2$fits$normal$coefficients))

  names(gf$fit$covariate$treat$coefficients) <-
    gsub("t0", "time", names(gf$fit$covariate$treat$coefficients))

  # outcome
  expect_equal(gf$fit$outcome$coefficients,
               gform_bin_eof$fits$outcome$coefficients)

  # covariates
  expect_equal(gf$fit$covariate$cov1$coefficients,
               gform_bin_eof$fits$cov1$coefficients)

  expect_equal(gf$fit$covariate$cov2$fits$zero$coefficients,
               gform_bin_eof$fits$cov2[[1]]$coefficients)

  expect_equal(gf$fit$covariate$cov2$fits$normal$coefficients,
               gform_bin_eof$fits$cov2[[2]]$coefficients)

  expect_equal(gf$fit$covariate$treat$coefficients,
               gform_bin_eof$fits$treat$coefficients)

})


test_that("scenario 5: natural course is the same", {
  nprisk_new <- s$np_means$obs_means
  nprisk_old <- gform_bin_eof$result

  expect_equal(
    nprisk_new$Py,
    nprisk_old$`NP mean`[nprisk_old$Interv. == 0]
  )
})

test_that("scenario 5: results are close", {
  result_new <- s$results
  result_old <- gform_bin_eof$result

  # risks
  expect_equal(
    result_new$means$estimate,
    result_old$`g-form mean`,
    tolerance = 0.02
  )

  result_new$ratios$estimate[is.na(result_new$ratios$estimate)] <- 1

  # ratios
  expect_equal(
    result_new$ratios$estimate,
    result_old$`Mean ratio`,
    tolerance = 0.02
  )

  result_new$diffs$estimate[is.na(result_new$diffs$estimate)] <- 0

  # differences
  expect_equal(
    result_new$diffs$estimate,
    result_old$`Mean difference`,
    tolerance = 0.02
  )
})


# scenario 6 --------------------------------------------------------------


