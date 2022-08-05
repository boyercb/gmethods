library(gfoRmula, exclude = c("gformula"), warn.conflicts = FALSE)
library(data.table)
library(readr)

# test visit process ------------------------------------------------------


# test 1 ------------------------------------------------------------------

# first verify that gmethods returns same results as gfoRmula without
# in HIV data without visit process

# gfoRmula package --------------------------------------------------------

dat_weeks <- read_csv("tests/testthat/practicum_data_weeks.csv", show_col_types = FALSE)

id <- "id"
time_name <- "week"
time_points <- 104 # two years
covnames <- c("visit", "pcp_new", "lwbc", "lcd4", "proph")
basecovs <- c("drug", "lcd4_0", "lwbc_0")
outcome_name <- "fail"
outcome_type <- "survival"
histories <- c(lagged, lagavg, cumavg)
histvars <- list(c("proph", "lcd4", "lwbc"), c("pcp_new"), c("pcp_new"))
covtypes <- c("binary", "binary", "normal", "normal", "binary")
covparams <- list(covmodels = c(
  visit ~ drug + lcd4_0 + lwbc_0 + lag1_proph + lag1_lcd4 + lag1_lwbc + lag_cumavg1_pcp_new + log(week + 1),
  pcp_new ~ drug + lcd4_0 + lwbc_0 + lag1_proph + lag1_lcd4 + lag1_lwbc + lag_cumavg1_pcp_new + log(week + 1),
  lwbc ~ drug + lcd4_0 + lwbc_0 + lag1_proph + lag1_lcd4 + lag1_lwbc + cumavg_pcp_new + log(week + 1),
  lcd4 ~ drug + lcd4_0 + lwbc_0 + lag1_proph + lag1_lcd4 + lwbc + cumavg_pcp_new + pcp_new + log(week + 1),
  proph ~ drug + lcd4_0 + lwbc_0 + lag1_proph + lcd4 + lwbc + cumavg_pcp_new + log(week + 1)
))

ymodel <- fail ~ drug*proph + lcd4_0 + lwbc_0 + lcd4 + lwbc + cumavg_pcp_new + log(week + 1)
intvars <- list(c("drug", "proph"), c("drug", "proph"), c("drug", "proph"), c("drug", "proph"))
int_times <- list(list(0:(time_points - 1), 0:(time_points - 1)), list(0:(time_points - 1), 0:(time_points - 1)),
                  list(0:(time_points - 1), 0:(time_points - 1)), list(0:(time_points - 1), 0:(time_points - 1)))
interventions <- list(
  list(c(static, rep(0, time_points)), c(static, rep(0, time_points))),
  list(c(static, rep(1, time_points)), c(static, rep(0, time_points))),
  list(c(static, rep(0, time_points)), c(static, rep(1, time_points))),
  list(c(static, rep(1, time_points)), c(static, rep(1, time_points)))
)

int_descript <- c(
  "drug = 0, proph = 0",
  "drug = 1, proph = 0",
  "drug = 0, proph = 1",
  "drug = 1, proph = 1"
)

dat_weeks <- dat_weeks[order(dat_weeks[[id]], dat_weeks[[time_name]]), ]
dat_weeks <- as.data.table(dat_weeks)

gform_res <- gfoRmula::gformula(
  obs_data = dat_weeks,
  id = id,
  time_name = time_name,
  time_points = time_points,
  covnames = covnames,
  basecovs = basecovs,
  outcome_name = outcome_name,
  outcome_type = outcome_type,
  covtypes = covtypes,
  covparams = covparams,
  ymodel = ymodel,
  intvars = intvars,
  interventions = interventions,
  int_descript = int_descript,
  int_times = int_times,
  ref_int = 1,
  histories = histories,
  histvars = histvars,
  sim_data_b = TRUE,
  nsimul = 5000,
  seed = runif(1, 0, .Machine$integer.max),
  model_fits = TRUE
)


# gmethods package --------------------------------------------------------

gf <- gformula(
  outcome_model = list(
    formula = fail ~ drug * proph + lcd4_0 + lwbc_0 + lcd4 + lwbc +
      cumavg_pcp_new + log(week + 1),
    link = "logit",
    family = "binomial"
  ),
  covariate_model = list(
    "pcp_new" = list(
      formula = pcp_new ~ drug + lcd4_0 + lwbc_0 + lag1_proph + lag1_lcd4 + lag1_lwbc +
        lag1_cumavg_pcp_new + log(week + 1),
      link = "logit",
      family = "binomial"
    ),
    "lwbc" = list(
      formula = lwbc ~ drug + lcd4_0 + lwbc_0 + lag1_proph + lag1_lcd4 +
        lag1_lwbc + cumavg_pcp_new + log(week + 1),
      link = "identity",
      family = "normal"
    ),
    "lcd4" = list(
      formula = lcd4 ~ drug + lcd4_0 + lwbc_0 + lag1_proph + lag1_lcd4 + lwbc +
        cumavg_pcp_new + pcp_new + log(week + 1),
      link = "identity",
      family = "normal"
    ),
    "proph" = list(
      formula = proph ~ drug + lcd4_0 + lwbc_0 + lag1_proph + lcd4 + lwbc +
        cumavg_pcp_new + log(week + 1),
      link = "logit",
      family = "binomial"
    )
  ),
  data = dat_weeks,
  survival = TRUE,
  id = 'id',
  time = 'week'
)

s <- simulate(
  gf,
  interventions = list(
    "drug = 0, proph = 0" = function(data, time) {
      set(data, j = "drug", value = 0)
      set(data, j = "proph", value = 0)
    },
    "drug = 1, proph = 0" = function(data, time) {
      set(data, j = "drug", value = 1)
      set(data, j = "proph", value = 0)
    },
    "drug = 0, proph = 1" = function(data, time) {
      set(data, j = "drug", value = 0)
      set(data, j = "proph", value = 1)
    },
    "drug = 1, proph = 1" = function(data, time) {
      set(data, j = "drug", value = 1)
      set(data, j = "proph", value = 1)
    }
  ),
  reference = "drug = 0, proph = 0",
  return_sims = TRUE,
  n_samples = 5000
)


# tests -------------------------------------------------------------------

# check that model coefficients are the same
test_that("no visit process: models fits are correct", {

  gfoRmula_fits <- gform_res$fits
  gmethods_fits <- gf$fit

  gfoRmula_fits <-
    lapply(gfoRmula_fits, function(x) {
      names(x$coefficients) <-
        gsub("lag_cumavg1_pcp_new",
             "lag1_cumavg_pcp_new",
             names(x$coefficients))
      return(x)
    })

  # outcome
  expect_equal(gmethods_fits$outcome$coefficients,
               gfoRmula_fits$fail$coefficients)

  # covariates
  expect_equal(gmethods_fits$covariate$pcp_new$coefficients,
               gfoRmula_fits$pcp_new$coefficients)

  expect_equal(gmethods_fits$covariate$lwbc$coefficients,
               gfoRmula_fits$lwbc$coefficients)

  expect_equal(gmethods_fits$covariate$lcd4$coefficients,
               gfoRmula_fits$lcd4$coefficients)

  expect_equal(gmethods_fits$covariate$proph$coefficients,
               gfoRmula_fits$proph$coefficients)

})

test_that("no visit process: natural course is the same", {
  nprisk_new <- s$np_means$obs_means
  nprisk_old <- gform_res$result

  expect_equal(
    nprisk_new$poprisk[nprisk_new$week == 103],
    nprisk_old$`NP Risk`[nprisk_old$k == 103 & nprisk_old$Interv. == 0]
  )
})

test_that("no visit process: results are close", {
  result_new <- s$results
  result_old <- gform_res$result[gform_res$result$k == 103, ]

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
    tolerance = 0.05
  )
})


# test 2 ------------------------------------------------------------------

# now verify that gmethods returns same results as gfoRmula in HIV
# data WITH visit process

# gfoRmula package --------------------------------------------------------

visit_process <- list(c('visit','pcp_new', Inf),
                      c('visit','lwbc', Inf) ,
                      c('visit','lcd4', Inf),
                      c('visit','proph', Inf))

gform_res <- gfoRmula::gformula(
  obs_data = dat_weeks,
  id = id,
  time_name = time_name,
  time_points = time_points,
  covnames = covnames,
  basecovs = basecovs,
  outcome_name = outcome_name,
  outcome_type = outcome_type,
  covtypes = covtypes,
  covparams = covparams,
  visitprocess = visit_process,
  ymodel = ymodel,
  intvars = intvars,
  interventions = interventions,
  int_descript = int_descript,
  int_times = int_times,
  ref_int = 1,
  histories = histories,
  histvars = histvars,
  sim_data_b = TRUE,
  seed = runif(1, 0, .Machine$integer.max),
  nsimul = 5000,
  model_fits = TRUE
)


# gmethods package --------------------------------------------------------

gf <- gformula(
  outcome_model = list(
    formula = fail ~ drug * proph + lcd4_0 + lwbc_0 + lcd4 + lwbc +
      cumavg_pcp_new + log(week + 1),
    link = "logit",
    family = "binomial"
  ),
  covariate_model = list(
    "pcp_new" = list(
      formula = pcp_new ~ drug + lcd4_0 + lwbc_0 + lag1_proph + lag1_lcd4 + lag1_lwbc +
        lag1_cumavg_pcp_new + log(week + 1),
      link = "logit",
      family = "binomial"
    ),
    "lwbc" = list(
      formula = lwbc ~ drug + lcd4_0 + lwbc_0 + lag1_proph + lag1_lcd4 +
        lag1_lwbc + cumavg_pcp_new + log(week + 1),
      link = "identity",
      family = "normal"
    ),
    "lcd4" = list(
      formula = lcd4 ~ drug + lcd4_0 + lwbc_0 + lag1_proph + lag1_lcd4 + lwbc +
        cumavg_pcp_new + pcp_new + log(week + 1),
      link = "identity",
      family = "normal"
    ),
    "proph" = list(
      formula = proph ~ drug + lcd4_0 + lwbc_0 + lag1_proph + lcd4 + lwbc +
        cumavg_pcp_new + log(week + 1),
      link = "logit",
      family = "binomial"
    )
  ),
  visit_model = list(
    formula = visit ~ drug + lcd4_0 + lwbc_0 + lag1_proph + lag1_lcd4 +
      lag1_lwbc + lag1_cumavg_pcp_new + log(week + 1),
    link = "logit",
    family = "binomial",
    covariates = c("pcp_new", "lwbc", "lcd4", "proph")
  ),
  data = dat_weeks,
  survival = TRUE,
  id = 'id',
  time = 'week'
)

s <- simulate(
  gf,
  interventions = list(
    "drug = 0, proph = 0" = function(data, time) {
      set(data, j = "drug", value = 0)
      set(data, j = "proph", value = 0)
    },
    "drug = 1, proph = 0" = function(data, time) {
      set(data, j = "drug", value = 1)
      set(data, j = "proph", value = 0)
    },
    "drug = 0, proph = 1" = function(data, time) {
      set(data, j = "drug", value = 0)
      set(data, j = "proph", value = 1)
    },
    "drug = 1, proph = 1" = function(data, time) {
      set(data, j = "drug", value = 1)
      set(data, j = "proph", value = 1)
    }
  ),
  reference = "drug = 0, proph = 0",
  bound_sims = TRUE,
  return_sims = TRUE,
  n_samples = 5000
)


# tests -------------------------------------------------------------------

# check that model coefficients are the same
test_that("visit process: models fits are correct", {

  gfoRmula_fits <- gform_res$fits
  gmethods_fits <- gf$fit

  gfoRmula_fits <-
    lapply(gfoRmula_fits, function(x) {
      names(x$coefficients) <-
        gsub("lag_cumavg1_pcp_new",
             "lag1_cumavg_pcp_new",
             names(x$coefficients))
      return(x)
    })

  # outcome
  expect_equal(gmethods_fits$outcome$coefficients,
               gfoRmula_fits$fail$coefficients)

  # covariates
  expect_equal(gmethods_fits$covariate$pcp_new$coefficients,
               gfoRmula_fits$pcp_new$coefficients)

  expect_equal(gmethods_fits$covariate$lwbc$coefficients,
               gfoRmula_fits$lwbc$coefficients)

  expect_equal(gmethods_fits$covariate$lcd4$coefficients,
               gfoRmula_fits$lcd4$coefficients)

  expect_equal(gmethods_fits$covariate$proph$coefficients,
               gfoRmula_fits$proph$coefficients)

  # visit model

  expect_equal(gmethods_fits$visit$coefficients,
               gfoRmula_fits$visit$coefficients)

})

test_that("visit process: natural course is the same", {
  nprisk_new <- s$np_means$obs_means
  nprisk_old <- gform_res$result

  expect_equal(
    nprisk_new$poprisk[nprisk_new$week == 103],
    nprisk_old$`NP Risk`[nprisk_old$k == 103 & nprisk_old$Interv. == 0]
  )
})

test_that("visit process: results are close", {
  result_new <- s$results
  result_old <- gform_res$result[gform_res$result$k == 103, ]

  # risks
  expect_equal(
    result_new$means$estimate,
    result_old$`g-form risk`,
    tolerance = 0.03
  )

  result_new$ratios$estimate[is.na(result_new$ratios$estimate)] <- 1

  # ratios
  expect_equal(
    result_new$ratios$estimate,
    result_old$`Risk ratio`,
    tolerance = 0.03
  )

  result_new$diffs$estimate[is.na(result_new$diffs$estimate)] <- 0

  # differences
  expect_equal(
    result_new$diffs$estimate,
    result_old$`Risk difference`,
    tolerance = 0.09
  )
})

