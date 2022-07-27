#' Define and fit (parametric) models for the g-formula
#'
#' @description
#' TBD
#'
#' @param outcome_model Model specification for outcome (Y).
#' @param covariate_model Model specifications for time-varying covariates,
#'   including potential treatment variables (L and A).
#' @param compevent_model Optional model specification for competing events (D).
#'   Only implemented for survival models.
#' @param censor_model Optional model specification for censoring (C). Only used
#'   to weight otherwise nonparametric estimates of the natural course such that
#'   the estimand is the same as those targeted by the (parametric) models for
#'   the outcome and covariates.
#' @param visit_model Optional model specification for visit process.
#' @param data A data.frame or matrix object containing the observed data.
#' @param id Name of unique id variable in `data`. Either a character vector or
#'   something coercible to one.
#' @param time Name of follow up time variable in `data`. Either a character
#'   vector or something coercible to one.
#' @param start_time Start of follow up time. Integer or scalar. Defaults to min
#'   value of `time`.
#' @param stop_time End of follow up time. Integer or scalar. Defaults to max
#'   value of `time`.
#' @param survival Outcome is a time-to-event/survival variable. Boolean.
#'   Default is FALSE.
#' @param censor_compevent Flag for whether competing events should be treated
#'   as censoring events (i.e. part of C) or modeled (i.e. part of L). Boolean.
#'   Only implemented for survival models. Default is TRUE.
#' @param max_visits If visit model specified, max number of missed visits
#'   before censoring (in model fit) or forcing visit (in simulation).
#' @param covs_baseline Names of time-fixed baseline covariates. Either a
#'   character vector or something coercible to one.
#' @param covs_tv Names of time-varying covariates. Either a
#'   character vector or something coercible to one.
#'
#' @return A gformula object with fitted models.
#' @import data.table
#' @export
#'
#' @examples
gformula <- function(outcome_model,
                     covariate_model,
                     compevent_model = NULL,
                     censor_model = NULL,
                     visit_model = NULL,
                     data,
                     id,
                     time,
                     start_time = 0,
                     stop_time = NULL,
                     survival = FALSE,
                     censor_compevent = TRUE,
                     max_visits = Inf,
                     covs_baseline = NULL,
                     covs_tv = NULL) {

  # check that data.table is installed
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package \"data.table\" must be installed to use this function.")
  }

  # If there's no time-varying covariates specified, use those listed in covariate model
  if (is.null(covs_tv)) {
    covs_tv <- names(covariate_model)
  } else {
    # If not present, add model covariates to list of time-varying covariates
    if (!all(covs_tv %in% names(covariate_model))) {
      covs_tv <- unique(c(covs_tv, names(covariate_model)))
    }
  }

  # Check that there's no overlap between baseline and time-varying covariates
  if (any(covs_tv %in% covs_baseline) & !(is.null(covs_tv) & is.null(covs_baseline))) {
    stop("Error: baseline and time-varying covariates must be distinct.")
  }

  # if not specified assume start time is first unique time
  if (is.null(start_time)) {
    start_time <- min(unique(data[[time]]))
  }

  # if not specified assume stop time is last unique time
  if (is.null(stop_time)) {
    stop_time <- max(unique(data[[time]]))
  }

  # add restrictions implied by visit process to covariate models
  if (!is.null(visit_model)) {
    f <- visit_model$formula
    visit <- all.vars(f)[attr(terms(f), "response")]

    # if user specifies covariates affected by visit process use those
    if (exists(visit_model, x = "covariates")) {
      visit_covs <- visit_model$covariates

      # otherwise assume all covariates are affected
    } else {
      visit_covs <- names(covariate_model)
    }

    covariate_model[names(covariate_model) %in% visit_covs] <-
      lapply(
        covariate_model[names(covariate_model) %in% visit_covs],
        function(object) {
          object$restrict <- list(
            subset = paste0(visit, " == 1")
          )
          return(object)
      })
  } else {
    visit_covs <- NULL
  }

  # get covariate restrictions
  restrictions <- lapply(covariate_model, function(object) {
    if (exists(object, x = "restrict")) {
      get(object, x = "restrict")
    } else {
      NULL
    }
  })



  # get list of all covariates that appear in model formulas (both dependent and
  # independent)
  covs <- lapply(
    c(covariate_model,
      list(outcome_model),
      if (!is.null(compevent_model)) list(compevent_model) else NULL,
      if (!is.null(censor_model)) list(censor_model) else NULL,
      if (!is.null(visit_model)) list(visit_model) else NULL
    ),
    function(object) {
    if (exists(object, x = "formula")) {
      all.vars(get(object, x = "formula"))
    } else {
      NULL
    }
  })
  covs <- unique(unlist(covs))

  if (!is.data.table(data)) {
    data <- as.data.table(data)
  } else {
    data <- copy(data)
  }

  # restrict to time range of interest
  rows <- which(data[[time]] >= start_time & data[[time]] <= stop_time)
  data <- data[rows, ]

  # create and document any history variables
  if (any(grepl("(^lag[0-9]+_)|(cumavg_)|(cumsum_)", covs))) {
    histories <- covs[grepl("(^lag[0-9]+_)|(cumavg_)|(cumsum_)", covs)]

    histories_to_add <- histories[!histories %in% colnames(data)]

    if (length(histories_to_add) > 0L) {
      make_histories(data, id, histories_to_add)
    }

    lagged_histories <- colnames(data)[grepl("^lag[0-9]+_", colnames(data))]
    cumulative_histories <- colnames(data)[grepl("(^cumavg_)|(^cumsum_)", colnames(data))]
  }

  # Fit outcome model
  outcome_fit <- fit_outcome_model(outcome_model, data, time, stop_time, survival)

  # Fit covariate models
  covariate_fit <- lapply(
    covariate_model,
    function(model) fit_covariate_model(model, data, time, start_time)
  )

  # Fit competing event model (if exists)
  if (!is.null(compevent_model)) {
    stopifnot(survival == TRUE)
    compevent_fit <- fit_compevent_model(compevent_model, data, time)
  } else {
    compevent_fit <- NULL
  }

  # Fit visit process model (if exists)
  if (!is.null(visit_model)) {
    visit_fit <- fit_covariate_model(visit_model, data, time, start_time)
  } else {
    visit_fit <- NULL
  }

  # Fit censoring model (if exists)
  if (!is.null(censor_model)) {
    censor_fit <- fit_compevent_model(censor_model, data, time)
  } else {
    censor_fit <- NULL
  }

  # Return g-formula definition object with fits
  ret <- list(
    outcome_model = outcome_model,
    covariate_model = covariate_model,
    compevent_model = compevent_model,
    censor_model = censor_model,
    visit_model = visit_model,
    fit = list(
      outcome = outcome_fit,
      covariate = covariate_fit,
      compevent = compevent_fit,
      censor = censor_fit,
      visit = visit_fit
    ),
    data = data,
    id = id,
    time = time,
    n_obs = length(unique(data[[id]])),
    start_time = start_time,
    stop_time = stop_time,
    survival = survival,
    censor_compevent = censor_compevent,
    max_visits = max_visits,
    visit_covs = visit_covs,
    covs_baseline = covs_baseline,
    covs_tv = covs_tv,
    lagged_histories = lagged_histories,
    cumulative_histories = cumulative_histories,
    restrictions = restrictions
  )

  class(ret) <- "gformula"

  return(ret)
}

#' Run Monte Carlo simulation to estimate population means under time-varying
#'   interventions using (parametric) g-formula
#'
#' @description
#' to approximate summation/integration over covariate history in g-formula
#'
#' @param object A g-formula object containing fitted models.
#' @param nsim Number of monte carlo samples to use to simulate histories.
#'   Integer. Defaults to number of unique IDs in either original data set or
#'   `newdata`.
#' @param seed Starting seed for random number generation. Integer.
#'   Default is runif(1, 0,  .Machine$integer.max).
#' @param newdata A data.frame or matrix object containing new data to use for
#'   simulation. Default is to use same dataset used to fit models. Useful for
#'   generalizing or transporting results to a target population under
#'   additional assumptions.
#' @param interventions A list of time-varying interventions
#' @param reference The name of the intervention to use as reference. Defaults
#'   to the natural course.
#' @param n_samples Number of monte carlo samples to use to simulate histories.
#'   Integer. Defaults to number of unique IDs in either original data set or
#'   `newdata`.
#' @param n_boots Number of bootstrap samples to draw. Integer. Defaults to `0`.
#' @param start_time Start of follow up time for simulation. Integer or scalar.
#'   Could be different than that specified during model fitting. Default is
#'   time specified in `gformula` fit.
#' @param stop_time End of follow up time for simulation. Integer or scalar.
#'   Could be different than that specified during model fitting. Default is
#'   time specified in `gformula` fit.
#' @param bound_sims Should covariate simulations be restricted to the range of
#'   the observed data? Boolean. Defaults to TRUE.
#' @param return_sims Flag for whether to return simulated covariate histories
#'   for each sample. Boolean. Defaults to FALSE.
#' @param natural_course Flag for whether to automatically add natural course to
#'   the list of interventions specified in `interventions`. Boolean. Defaults
#'   to TRUE.
#' @param conf_level A scalar specifying the confidence level for bootstrapped
#'   intervals in output.
#' @param conf_type A string specifying the type of bootstrap confidence
#'   interval requested. One of "norm" or "perc".
#'
#' @return a gformula.simulation object with
#' @import data.table
#' @export
#'
#' @examples
simulate.gformula <- function(object,
                              nsim = NULL,
                              seed = runif(1, 0, .Machine$integer.max),
                              newdata = NULL,
                              interventions = NULL,
                              reference = "Natural course",
                              n_samples = NULL,
                              n_boots = NULL,
                              start_time = NULL,
                              stop_time = NULL,
                              bound_sims = TRUE,
                              return_sims = FALSE,
                              natural_course = TRUE,
                              conf_level = 0.95,
                              conf_type = "norm",
                              ...) {
  set.seed(seed)

  # add natural course to intervention list
  if (is.null(interventions)) {
    interventions <- list(
      "Natural course" = "natural course"
    )
  } else {
    if (!is.list(interventions)) {
      interventions <- list(interventions)
    } else if (!"natural course" %in% interventions & natural_course == TRUE) {
      interventions <- c(list("Natural course" = "natural course"), interventions)
    }
  }

  # update start and stop time
  if (is.null(start_time)) {
    start_time <- object$start_time
  }

  if (is.null(stop_time)) {
    stop_time <- object$stop_time
  }

  if (!reference %in% names(interventions)) {
    stop("Reference intervention must correspond to named intervention.")
  }

  res <- run_gformula(
    gformula = object,
    newdata = newdata,
    interventions = interventions,
    n_sims = 1,
    n_samples = n_samples,
    n_boots = n_boots,
    start_time = start_time,
    stop_time = stop_time,
    bound_sims = bound_sims,
    return_sims = return_sims,
    last_only = TRUE,
    natural_course = natural_course,
    conf_level = conf_level,
    conf_type = conf_type
  )

  sims <- res$sims
  cov_means <- res$covs
  out_means <- res$means
  boots <- res$boots

  t0 <- res$t0
  t <- res$t

  target <- switch(object$fit$outcome$type,
    "survival" = "poprisk",
    "continuous" = "Ey",
    "binomial" = "Py"
  )

  # if there are bootstrap estimates
  if (!is.null(n_boots)) {
    if (object$fit$outcome$type == "survival") {
      # take only last time point for summary
      means <- lapply(boots$means, function(x) x[nrow(x), ])
      t0 <- lapply(t0$means, function(x) x[length(x), ][[target]])
      t <- lapply(t$means, function(x) x[, ncol(x)])
    } else {
      means <- boots$means
      t0 <- lapply(t0$means, function(x) x[[target]])
      t <- t$means
    }

    out_means <- boots$means
    cov_means <- boots$covs

    means <- do.call("rbind", means)
    means[[object$time]] <- NULL
  } else {
    if (object$fit$outcome$type == "survival") {
      means <- lapply(out_means, function(x) x[nrow(x), ])
    } else {
      means <- out_means
    }

    means <- do.call("rbind", means)

    # add NAs for confidence limits
    means <- cbind(means, NA, NA)
    means[[object$time]] <- NULL
    colnames(means) <- c("estimate", "conf.low", "conf.high")
  }

  # Hold special place for natural course in intervention list (0)
  if ("Natural course" %in% names(interventions)) {
    means <- cbind("intervention" = 1:nrow(means) - 1, means)
  } else {
    means <- cbind("intervention" = 1:nrow(means), means)
  }

  means <- as.data.frame(means)

  ref <- which(names(interventions) %in% reference)
  nonref <- which(!names(interventions) %in% reference)

  diffs <- means
  if (!is.null(n_boots)) {
    diffs_ci <-
      lapply(nonref, function(i) {
        get_bootstrap_ci(
          t = t[[i]] - t[[ref]],
          t0 = t0[[i]] - t0[[ref]],
          conf_level = conf_level,
          conf_type = conf_type
        )
      })
    diffs[-ref, 2:4] <- do.call("rbind", diffs_ci)
  } else {
    diffs[-ref, 2] <- diffs[-ref, 2] - diffs[ref, 2]
  }
  diffs[ref, 2:4] <- NA

  ratios <- means
  if (!is.null(n_boots)) {
    ratios_ci <-
      lapply(nonref, function(i) {
        get_bootstrap_ci(
          t = t[[i]] / t[[ref]],
          t0 = t0[[i]] / t0[[ref]],
          conf_level = conf_level,
          conf_type = conf_type
        )
      })
    ratios[-ref, 2:4] <- do.call("rbind", ratios_ci)
  } else {
    ratios[-ref, 2] <- ratios[-ref, 2] / ratios[ref, 2]
  }
  ratios[ref, 2:4] <- NA

  np_means <- calc_np_means(
    outcome_fit = object$fit$outcome,
    covariate_fit = object$fit$covariate,
    compevent_fit = object$fit$compevent,
    censor_fit = object$fit$censor,
    visit_fit = object$fit$visit,
    data = if (!is.null(newdata)) newdata else object$data,
    id = object$id,
    time = object$time,
    censor_compevent = object$censor_compevent,
    stop_time = stop_time
  )

  ret <- list(
    results = list(
      means = means,
      ratios = ratios,
      diffs = diffs
    ),
    sims = sims,
    boots = list(
      t0 = res$t0,
      t = res$t
    ),
    out_means = lapply(out_means, as.data.frame),
    cov_means = lapply(cov_means, as.data.frame),
    np_means = np_means,
    interventions = interventions,
    type = object$fit$outcome$type,
    reference = reference,
    id = object$id,
    time = object$time,
    start_time = start_time,
    stop_time = stop_time,
    n_obs = object$n_obs,
    n_samples = n_samples,
    n_boots = n_boots,
    conf_level = conf_level
  )


  class(ret) <- "gformula.simulation"

  return(ret)
}

#' Predict conditional expectation for specified covariate profile under
#'   time-varying interventions using the g-formula
#'
#' @description
#' TBD.
#'
#' @param object A g-formula object containing fitted models.
#' @param newdata A data.frame or matrix object containing new data to use for
#'   simulation. Default is to use same dataset used to fit models. Useful for
#'   making out-of-sample predictions.
#' @param interventions A list of time-varying interventions
#' @param n_sims Number of monte carlo simulations PER individual. For example,
#'   for a dataset with 2000 unique individuals (as determined by identifier
#'   specified for `id`) an `n_sims` of 100 will perform 100 simulatios for
#'   each observation to generate an "individual" prediction.
#'   Integer. Defaults to 100.
#' @param n_boots Number of bootstrap samples to draw. Integer. Defaults to `0`.
#' @param start_time Start of follow up time for simulation. Integer or scalar.
#'   Could be different than that specified during model fitting. Default is
#'   time specified in `gformula` fit.
#' @param stop_time End of follow up time for simulation. Integer or scalar.
#'   Could be different than that specified during model fitting. Default is
#'   time specified in `gformula` fit.
#' @param seed Starting seed for random number generation.
#' @param bound_sims Should covariate simulations be restricted to the range of
#'   the observed data? Boolean. Defaults to FALSE.
#' @param return_sims Flag for whether to return simulated covariate histories
#'   for each sample. Boolean. Defaults to FALSE.
#' @param natural_course
#' @param last_only Flag for whether prediction should be for final time point
#'   only (Default) or should predictions at intermediate time points be
#'   included? Boolean. Most relevant for survival models.
#' @param predict_covs Flag for whether covariate predictions should be included
#'   in output. Boolean. Default is False.
#' @param conf_level A scalar specifying the confidence level for bootstrapped
#'   intervals in output.
#' @param conf_type A string specifying the type of bootstrap confidence
#'   interval requested. One of "norm" or "perc".
#'
#' @return
#' @export
#'
#' @examples
predict.gformula <- function(object,
                             newdata = NULL,
                             interventions = NULL,
                             n_sims = 100,
                             n_boots = NULL,
                             start_time = NULL,
                             stop_time = NULL,
                             seed = runif(1, 0, .Machine$integer.max),
                             bound_sims = TRUE,
                             return_sims = FALSE,
                             natural_course = TRUE,
                             last_only = TRUE,
                             predict_covs = FALSE,
                             conf_level = 0.95,
                             conf_type = "norm",
                             ...) {
  set.seed(seed)

  # update start and stop time
  if (is.null(start_time)) {
    start_time <- object$start_time
  }

  if (is.null(stop_time)) {
    stop_time <- object$stop_time
  }

  # add natural course to intervention list
  if (is.null(interventions)) {
    interventions <- list(
      "Natural course" = "natural course"
    )
  } else {
    if (!is.list(interventions)) {
      interventions <- list(interventions)
    } else if (!"natural course" %in% interventions &
      natural_course == TRUE) {
      interventions <-
        c(list("Natural course" = "natural course"), interventions)
    }
  }

  res <- run_gformula(
    gformula = object,
    newdata = newdata,
    interventions = interventions,
    n_sims = n_sims,
    n_samples = NULL,
    n_boots = n_boots,
    start_time = start_time,
    stop_time = stop_time,
    bound_sims = bound_sims,
    return_sims = return_sims,
    last_only = last_only,
    natural_course = FALSE,
    conf_level = conf_level,
    conf_type = conf_type
  )

  if (!is.null(n_boots)) {
    preds <- res$boots$means
    covs <- res$boots$covs
  } else {
    preds <- res$means
    covs <- res$covs
  }

  if (predict_covs) {
    # join on id, time
    preds <- lapply(seq_along(covs), function(i) {
      preds[[i]][covs[[i]], on = c(object$id, object$time)]
    })

    names(preds) <- names(covs)
  }

  # convert to data.frame
  preds <- lapply(preds, as.data.frame)

  if (length(preds) == 1) {
    preds <- preds[[1]]
  } else {
    class(preds) <- "gformula.prediction"

  }

  return(preds)
}



# s3 helpers --------------------------------------------------------------


#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
# print.gformula <- function(x) {
#
# }


#' @export
print.gformula.simulation <- function(x, ...) {
  if (!inherits(x, "gformula.simulation")) {
    stop("Argument 'x' must be an object of class \"gformula.simulation\".")
  }

  header <- get_header(
    ints = names(x$interventions),
    ref = which(names(x$interventions) %in% x$reference)
  )

  footer = get_footer(
    size = x$n_obs,
    n_samples = if (!is.null(x$n_samples)) x$n_samples else x$n_obs,
    n_boots = if (!is.null(x$n_boots)) x$n_boots else 0
  )

  cat(header, "\n\n")

  if (!is.null(x$n_boots)) {
    colnames(x$results$means) <-
      c("Intervention",
        switch(
          x$type,
          "survival" = "Risk",
          "continuous" = "Mean",
          "binomial" = "Risk"
        ),
        paste0(round(100 * x$conf_level), "% CI lower"),
        paste0(round(100 * x$conf_level), "% CI upper"))

    colnames(x$results$diffs) <-
      c("Intervention",
        switch(
          x$type,
          "survival" = "RD",
          "continuous" = "Mean Diff",
          "binomial" = "RD"
        ),
        paste0(round(100 * x$conf_level), "% CI lower"),
        paste0(round(100 * x$conf_level), "% CI upper"))

    colnames(x$results$ratios) <-
      c("Intervention",
        switch(
          x$type,
          "survival" = "RR",
          "continuous" = "Mean Ratio",
          "binomial" = "RR"
        ),
        paste0(round(100 * x$conf_level), "% CI lower"),
        paste0(round(100 * x$conf_level), "% CI upper"))

    cat(" Population Means:\n")

    print(x$results$means, row.names = FALSE, col.names='top', digits = 4)

    cat("\n Ratios:\n")

    ratios <- format.data.frame(x$results$ratios, digits = 4, trim = TRUE)
    ratios[ratios == "NA"] <- "ref"
    print(ratios, row.names = FALSE, col.names='top')

    cat("\n Differences:\n")

    diffs <- format.data.frame(x$results$diffs, digits = 4, trim = TRUE)
    diffs[diffs == "NA"] <- "ref"
    print(diffs, row.names = FALSE, col.names='top')

    cat("\n\n Natural Course Benchmarks:\n")
    if (x$type == "survival") {
      nc <- merge(
        x$np_means$obs_means,
        x$out_means[["Natural course"]],
        by = x$time
      )
    } else {
      nc <- cbind(
        x$stop_time,
        x$np_means$obs_means,
        x$out_means[["Natural course"]][, -1]
      )
    }
    nc <- as.data.frame(nc)
    colnames(nc) <- c("Time",
                      switch(
                        x$type,
                        "survival" = "NP Risk",
                        "continuous" = "NP Mean",
                        "binomial" = "NP Risk"
                      ),
                      switch(
                        x$type,
                        "survival" = "G-formula Risk",
                        "continuous" = "G-formula Mean",
                        "binomial" = "G-formula Risk"
                      ),
                      paste0(round(100 * x$conf_level), "% CI lower"),
                      paste0(round(100 * x$conf_level), "% CI upper"))
    nc <- format.data.frame(nc, digits = 4, trim = TRUE)
    print(nc, row.names = FALSE, col.names='top')
  } else {
    result <- cbind(
      x$results$means[, 1:2],
      x$results$diffs[, 2],
      x$results$ratios[, 2]
    )

    colnames(result) <-
      c("Intervention",
        switch(
          x$type,
          "survival" = "Risk",
          "continuous" = "Mean",
          "binomial" = "Risk"
        ),
        switch(
          x$type,
          "survival" = "RD",
          "continuous" = "Mean Diff",
          "binomial" = "RD"
        ),
        switch(
          x$type,
          "survival" = "RR",
          "continuous" = "Mean Ratio",
          "binomial" = "RR"
        ))

    cat(" Results:\n")

    result <- format.data.frame(result, digits = 4, trim = TRUE)
    result[result == "NA"] <- "ref"
    print(result, row.names = FALSE, col.names='top')

    cat("\n\n Natural Course Benchmarks:\n")
    if (x$type == "survival") {
      nc <- merge(
        x$np_means$obs_means,
        x$out_means[["Natural course"]],
        by = x$time
      )
    } else {
      nc <- cbind(
        x$stop_time,
        x$np_means$obs_means,
        x$out_means[["Natural course"]][, -1]
      )
    }
    nc <- as.data.frame(nc)
    colnames(nc) <- c("Time",
                      switch(
                        x$type,
                        "survival" = "NP Risk",
                        "continuous" = "NP Mean",
                        "binomial" = "NP Risk"
                      ),
                      switch(
                        x$type,
                        "survival" = "G-formula Risk",
                        "continuous" = "G-formula Mean",
                        "binomial" = "G-formula Risk"
                      ))
    nc <- format.data.frame(nc, digits = 4, trim = TRUE)
    print(nc, row.names = FALSE, col.names='top')

  }

  cat(paste0("\n\nStart time: ", x$start_time, ", Stop time: ", x$stop_time))
  cat("\n", footer)

}


#' @export
summary.gformula <- function(object, ...) {
  if (!inherits(object, "gformula")){
    stop("Argument 'object' must be an object of class \"gformula\".")
  }
  class(object) <- c("summary.gformula", class(object))
  return (object)
}

#' @export
summary.gformula.simulation <- function(object, ...) {
  if (!inherits(object, "gformula.simulation")){
    stop("Argument 'object' must be an object of class \"gformula.simulation\".")
  }
  class(object) <- c("summary.gformula.simulation", class(object))
  return (object)
}

#' @export
plot.gformula <- function(x, ...) {

}

#' @export
plot.gformula.simulation <- function(x, ...) {

}


# internal helpers --------------------------------------------------------

#' Title
#'
#' @param gformula
#' @param newdata
#' @param interventions
#' @param n_sims
#' @param n_samples
#' @param n_boots
#' @param start_time
#' @param stop_time
#' @param bound_sims
#' @param return_sims
#' @param last_only
#' @param natural_course
#' @param conf_level
#' @param conf_type
#'
#' @return
#'
#' @examples
run_gformula <- function (
    gformula,
    newdata = NULL,
    interventions = NULL,
    n_sims = 1,
    n_samples = NULL,
    n_boots = NULL,
    start_time = NULL,
    stop_time = NULL,
    bound_sims = FALSE,
    return_sims = FALSE,
    last_only = TRUE,
    natural_course = TRUE,
    conf_level = 0.95,
    conf_type = "norm"
) {

  # copy variables from g-formula object
  id <- gformula$id
  time <- gformula$time
  covs_baseline <- gformula$covs_baseline
  covs_tv <- gformula$covs_tv
  lagged_histories <- gformula$lagged_histories
  cumulative_histories <- gformula$cumulative_histories
  restrictions <- gformula$restrictions

  if (is.null(newdata)) {
    data <- gformula$data
  } else {
    data <- newdata

    # maybe add some error checking here!
  }

  covariate_fit <- gformula$fit$covariate
  outcome_fit <- gformula$fit$outcome
  compevent_fit <- gformula$fit$compevent
  censor_fit <- gformula$fit$censor
  visit_fit <- gformula$fit$visit

  # local variables
  uids <- unique(data[[id]])
  utimes <- unique(data[[time]])

  # collect variable names
  covs <- names(covariate_fit)
  outcome <- outcome_fit$outcome

  if (!is.null(compevent_fit)) {
    compevent <- compevent_fit$outcome
  } else {
    compevent <- NULL
  }

  if (!is.list(interventions)) {
    interventions <- list(interventions)
  }

  types <- lapply(covariate_fit, get, x = "type")

  if (bound_sims) {
    # Determine ranges of observed covariates and outcome
    covariate_range <- lapply(seq_along(covs), function(i) {
      if (types[i] %in% c('normal',
                          'bounded normal',
                          'truncated normal',
                          'zero-inflated normal',
                          'poisson',
                          'zero-inflated poisson',
                          'zero-inflated negbin',
                          'poisson hurdle',
                          'negbin hurdle'
      )) {
        range(gformula$data[[covs[i]]], na.rm = TRUE)
      } else {
        NULL
      }
    })
    outcome_range <- range(gformula$data[[outcome]], na.rm = TRUE)

    if (!is.null(compevent_fit)) {
      compevent_range <- range(gformula$data[[compevent]])
    } else {
      compevent_range <- NULL
    }
  } else {
    covariate_range <- NULL
    outcome_range <- NULL
    compevent_range <- NULL
  }

  visit_covs <- gformula$visit_covs

  # initialize data for data.table
  rows <- which(data[[time]] == start_time)
  dt <- data[rows, ]

  # add hidden ordering variable
  obs <- nrow(dt)
  dt$.order <- 1:obs

  # create data.table
  setDT(dt)

  # If the number of desired simulation units is greater than the number of
  # unique units sample with replacement
  if (!is.null(n_samples)) {
    if (n_samples != obs) {
      dt <- dt[sample(.N, n_samples, replace = TRUE)]
    }
  }

  dt[, .sid := 1:nrow(.SD)]

  # If number of simulations PER unit is greater than 1
  if (n_sims > 1) {
    prediction <- TRUE
    # make copies
    dt <- dt[rep(1:.N, each = n_sims)]
    dt[, .sim := 1:.N, by = .sid]
  } else {
    prediction <- FALSE
    dt[, .sim := 1]
  }

  # make sure data is sorted by id and sim
  setorderv(dt, c(id, '.sid', '.sim'))

  sim_results <-
    lapply(
      interventions,
      function (x) {
        # simulate interventions
        simulate_intervention(
          dt,
          x,
          id,
          time,
          start_time,
          stop_time,
          outcome_fit,
          covariate_fit,
          compevent_fit,
          censor_fit,
          visit_fit,
          visit_covs,
          lagged_histories,
          cumulative_histories,
          covariate_range,
          outcome_range,
          compevent_range,
          restrictions,
          bound_sims,
          return_sims,
          last_only,
          prediction
        )
      })

  # extract results
  sims <- lapply(sim_results, get, x = "sims")
  means <- lapply(sim_results, get, x = "means")
  covs <- lapply(sim_results, get, x = "covs")

  names(sims) <- names(interventions)
  names(means) <- names(interventions)
  names(covs) <- names(interventions)

  # Use nonparametric bootstrap to estimate uncertainty
  if (!is.null(n_boots)) {

    # create splits based on time blocks
    splits <- split(1:nrow(gformula$data), gformula$data[[id]])

    target <- switch(
      gformula$fit$outcome$type,
      "survival" = "poprisk",
      "continuous" = "Ey",
      "binomial" = "Py"
    )

    key <- construct_key(means, covs, target, interventions, prediction)

    # calculate bootstrap statistics
    bs <- boot::boot(
      data = uids,
      statistic = bootstrap_simulations,
      R = n_boots,
      gformula = gformula,
      splits = splits,
      newdata = newdata,
      interventions = interventions,
      n_sims = n_sims,
      n_samples = n_samples,
      start_time = start_time,
      stop_time = stop_time,
      bound_sims = bound_sims,
      prediction = prediction,
      last_only = last_only
    )

    ests <- get_bootstrap_ci(bs$t, bs$t0, conf_level, conf_type)

    boots <- list(
      covs = use_covs_key(key, covs, ests, prediction),
      means = use_means_key(key, means, ests, prediction)
    )

    t0 <- list(
      covs = use_covs_key(key, covs, bs$t0, prediction),
      means = use_means_key(key, means, bs$t0, prediction)
    )

    t <- list(
      covs = use_covs_key_b(key, covs, bs$t, prediction),
      means = use_means_key_b(key, means, bs$t, prediction)
    )

  } else {
    boots <- NULL
    t0 <- NULL
    t <- NULL
  }

  ret <- list(
    sims = sims,
    means = means,
    covs = covs,
    boots = boots,
    t0 = t0,
    t = t
  )

  return(ret)
}

#' Title
#'
#' @param dt
#' @param intervention
#' @param id
#' @param time
#' @param start_time
#' @param stop_time
#' @param outcome_fit
#' @param covariate_fit
#' @param compevent_fit
#' @param censor_fit
#' @param visit_fit
#' @param visit_covs
#' @param lagged_histories
#' @param cumulative_histories
#' @param covariate_range
#' @param outcome_range
#' @param compevent_range
#' @param restrictions
#' @param bound_sims
#' @param return_sims
#' @param last_only
#' @param prediction
#'
#' @return
#'
#' @examples
simulate_intervention <-
  function(dt,
           intervention,
           id,
           time,
           start_time,
           stop_time,
           outcome_fit,
           covariate_fit,
           compevent_fit,
           censor_fit,
           visit_fit,
           visit_covs,
           lagged_histories,
           cumulative_histories,
           covariate_range,
           outcome_range,
           compevent_range,
           restrictions,
           bound_sims,
           return_sims,
           last_only,
           prediction) {

    nsims <- nrow(dt)
    sims <- list(dt)

    covs <- names(covariate_fit)
    rmses <- lapply(covariate_fit, get, x = "rmse")
    types <- lapply(covariate_fit, get, x = "type")

    if (!is.null(compevent_fit)) {
      compevent <- compevent_fit$outcome
    } else {
      compevent <- NULL
    }

    if (!is.null(visit_fit)) {
      visit <- visit_fit$outcome
    } else {
      visit <- NULL
    }

    for (k in start_time:stop_time) {
      i <- which(start_time:stop_time == k)

      if (k  == start_time) {

        # initialize current sim
        sim <- copy(sims[[i]])

      } else {
        # set initial values of sim to previous
        sim <- copy(sims[[i - 1]])

        # update time
        set(sim, j = time, value = rep(k, nsims))

        # update lags in covariate histories
        update_lagged_histories(sim, lagged_histories)

        # if visit model specified determine whether visit occurs
        if (!is.null(visit_fit)) {
          set(
            sim,
            j = visit,
            value = draw_binomial(
              N = nsims,
              size = 1,
              p.fit = visit_fit,
              data = sim
            )
          )
        }

        # draw new covariate values
        for (j in seq_along(covs)) {

          # visit process variables, default is to draw values for everyone
          visit_rows <- 1:nrow(sim)
          ndraws <- nsims

          # if covariate is part of visit process draw values only at visits
          if (!is.null(visit_fit)) {
            if (!covs[[j]] %in% visit_covs) {
              visit_rows <- which(sim[[visit]] == 1)
              ndraws <- length(visit_rows)
            }
          }

          if (types[[j]] == "binomial") {
            set(
              sim,
              i = visit_rows,
              j = covs[[j]],
              value = draw_binomial(
                N = ndraws,
                size = 1,
                p.fit = covariate_fit[[j]],
                data = sim[visit_rows, ]
              )
            )
          } else if (types[[j]] == 'normal') {
            set(
              sim,
              i = visit_rows,
              j = covs[[j]],
              value = draw_normal(
                N = ndraws,
                mu.fit = covariate_fit[[j]],
                sd.hat = rmses[[j]],
                data = sim[visit_rows, ]
              ))

          } else if (types[[j]] == 'poisson') {
            set(
              sim,
              i = visit_rows,
              j = covs[[j]],
              value = draw_poisson(
                N = ndraws,
                lambda.fit = covariate_fit[[j]],
                data = sim[visit_rows, ]
              ))
          } else if (types[[j]] == 'categorical') {
            set(
              sim,
              i = visit_rows,
              j = covs[[j]],
              value = draw_multinomial(
                p.fit = covariate_fit[[j]],
                data = sim[visit_rows, ]
              ))
          } else if (types[[j]] == 'zero-inflated normal') {
            set(
              sim,
              i = visit_rows,
              j = covs[[j]],
              value = draw_zero_inflated_normal(
                N = ndraws,
                size = 1,
                p.fit = covariate_fit[[j]]$fits$zero,
                mu.fit = covariate_fit[[j]]$fits$normal,
                sd.hat = rmses[[j]],
                data = sim[visit_rows, ]
              ))
          } else if (types[[j]] == 'zero-inflated poisson') {
            set(
              sim,
              i = visit_rows,
              j = covs[[j]],
              value = draw_zero_inflated_poisson(
                N = ndraws,
                size = 1,
                fit = covariate_fit[[j]],
                data = sim[visit_rows, ]
              ))
          } else if (types[[j]] == 'zero-inflated negbin') {
            set(
              sim,
              i = visit_rows,
              j = covs[[j]],
              value = draw_zero_inflated_negbin(
                N = ndraws,
                size = 1,
                fit = covariate_fit[[j]],
                data = sim[visit_rows, ]
              ))
          } else if (types[[j]] == 'poisson hurdle') {
            set(
              sim,
              i = visit_rows,
              j = covs[[j]],
              value = draw_poisson_hurdle(
                N = ndraws,
                size = 1,
                fit = covariate_fit[[j]],
                data = sim[visit_rows, ]
              ))
          } else if (types[[j]] == 'negbin hurdle') {
            set(
              sim,
              i = visit_rows,
              j = covs[[j]],
              value = draw_negbin_hurdle(
                N = ndraws,
                size = 1,
                fit = covariate_fit[[j]],
                data = sim[visit_rows, ]
              ))
          } else if (types[[j]] == 'bounded normal') {
            set(
              sim,
              i = visit_rows,
              j = covs[[j]],
              value = draw_bounded_normal(
                N = ndraws,
                mu.fit = covariate_fit[[j]],
                sd.hat = rmses[[j]],
                range = covariate_range[[j]],
                data = sim[visit_rows, ]
              ))
          } else if (types[[j]] == 'truncated normal') {
            set(
              sim,
              i = visit_rows,
              j = covs[[j]],
              value = draw_truncated_normal(
                N = ndraws,
                mu.fit = covariate_fit[[j]],
                sd.hat = rmses[[j]],
                direction = covariate_fit[[j]]$direction,
                point = covariate_fit[[j]]$point,
                data = sim[visit_rows, ]
              ))
          } else if (types[[j]] == 'binomial grf') {
            set(
              sim,
              i = visit_rows,
              j = covs[[j]],
              value = draw_grf_binomial(
                N = ndraws,
                size = 1,
                p.fit = covariate_fit[[j]],
                data = sim[visit_rows, covariate_fit[[j]]$covs]
              ))
          } else if (types[[j]] == 'normal grf') {
            set(
              sim,
              i = visit_rows,
              j = covs[[j]],
              value = draw_grf_normal(
                N = ndraws,
                mu.fit = covariate_fit[[j]],
                sd.hat = rmses[[j]],
                data = sim[visit_rows, covariate_fit[[j]]$covs]
              ))
          } else if (types[[j]] == 'custom') {
            set(
              sim,
              i = visit_rows,
              j = covs[[j]],
              value = draw_custom(
                N = ndraws,
                mu.fit = covariate_fit[[j]],
                sd.hat = rmses[[j]],
                data = sim[visit_rows, ],
                pred.fun = custom.pred.fun
              ))
          }

          if (bound_sims) {
            # Set simulated covariate values outside the observed range to the observed min / max
            if (types[[j]] %in% c('normal', 'bounded normal', 'truncated normal')) {
              if (min(sim[[covs[j]]]) < covariate_range[[j]][1]) {
                sim[sim[[covs[j]]] < covariate_range[[j]][1],
                    (covs[j]) := covariate_range[[j]][1]]
              }
              if (max(sim[[covs[j]]]) > covariate_range[[j]][2]) {
                sim[sim[[covs[j]]] > covariate_range[[j]][2],
                    (covs[j]) := covariate_range[[j]][2]]
              }

            } else if (types[[j]] %in% c(
              'zero-inflated normal',
              'zero-inflated poisson',
              'zero-inflated negbin',
              'poisson hurdle',
              'negbin hurdle'
            )) {
              if (min(sim[[covs[j]]]) < covariate_range[[j]][1]) {
                sim[sim[[covs[j]]] < covariate_range[[j]][1] &
                      sim[[covs[j]]] > 0, (covs[j]) := covariate_range[[j]][1]]
              }
              if (max(sim[[covs[j]]]) > covariate_range[[j]][2]) {
                sim[sim[[covs[j]]] > covariate_range[[j]][2],
                    (covs[j]) := covariate_range[[j]][2]]
              }
            }
          }

          # restrictions
          if (!is.null(restrictions[[j]])) {
            if (exists(restrictions[[j]], x = "otherwise")) {
              rows <- with(sim, !eval(parse(text = restrictions[[j]]$subset)))
              sim[[covs[j]]][rows] <- restrictions[[j]]$otherwise
            }
          }

          # if covariate is part of cumulative history update
          if (covs[[j]] %in% gsub("^cumsum_", "", cumulative_histories)) {
            update_cumulative_histories(sim, sims, paste0("cumsum_", covs[[j]]))
          } else if (covs[[j]] %in% gsub("^cumavg_", "", cumulative_histories)) {
            update_cumulative_histories(sim, sims, paste0("cumavg_", covs[[j]]))
          }
        }
      }

      # intervene on treatment variables
      if (is.function(intervention)) {
        if (!prediction) {
          intervention(sim, k)
        } else {
          # predictions are assumed to be conditional on values at start time
          # TODO: think about whether this needs to be imposed or just explained
          if (k > start_time) {
            intervention(sim, k)
          }
        }

      } else if (intervention != "natural course") {
        stop(paste("Error don't know what to do with intervention:", intervention))
      }

      # Predict outcome value at time t using parametric models
      if (outcome_fit$type == 'survival') {
        set(
          sim,
          j = 'Py',
          value = stats::predict(outcome_fit, type = 'response', newdata = sim)
        )

        # Predict competing event probabilities
        if (!is.null(compevent)) {
          set(
            sim,
            j = 'Pd',
            value = stats::predict(compevent_fit, type = 'response', newdata = sim)
          )

          # Simulate competing event variable
          set(sim, j = 'D', value = stats::rbinom(nsims, 1, sim$Pd))

          # Set simulated compevent values outside the observed range to the observed min / max
          if (bound_sims) {
            if (length(sim[sim$D < compevent_range[1]]$D) != 0){
              sim[sim$D < compevent_range[1], 'D' := compevent_range[1]]
            }
            if (length(sim[sim$D > compevent_range[2]]$D) != 0){
              sim[sim$D > compevent_range[2], 'D' := compevent_range[2]]
            }
          }
        } else {
          set(sim, j = 'D', value = 0)
          set(sim, j = 'Pd', value = 0)
        }
      } else if (outcome_fit$type == 'continuous') {
        if (k < stop_time) {
          set(sim, j = 'Ey', value = as.double(NA))
        } else if (k == stop_time) {
          set(
            sim,
            j = 'Ey',
            value = stats::predict(outcome_fit, type = 'response', newdata = sim)
          )
        }
      } else if (outcome_fit$type == 'binomial') {
        if (k < stop_time) {
          set(sim, j = 'Py', value = as.double(NA))
        } else if (k == stop_time) {
          set(
            sim,
            j = 'Py',
            value = stats::predict(outcome_fit, type = 'response', newdata = sim)
          )
        }
      }

      # Simulate outcome variable
      if (outcome_fit$type == 'survival') {
        set(sim, j = 'Y', value = stats::rbinom(nsims, 1, sim$Py))

        # If competing event occurs, outcome cannot also occur because
        # both presumably lead to death
        sim[sim$D == 1, 'Y' := NA]
      }

      # Set simulated outcome values outside the observed range to the observed min / max
      if (bound_sims) {
        if (length(sim[sim$Y < outcome_range[1]]$Y) != 0){
          sim[sim$Y < outcome_range[1], 'Y' := outcome_range[1]]
        }
        if (length(sim[sim$Y > outcome_range[2]]$Y) != 0){
          sim[sim$Y > outcome_range[2], 'Y' := outcome_range[2]]
        }
      }

      if (outcome_fit$type == 'survival') {
        if (k == start_time) {
          # Calculate probability of Y rather than D at time t
          set(sim, j = 'prodp1', value = sim$Py * (1 - sim$Pd))

          # Calculate probability of survival or D = 1
          set(sim, j = 'prodp0', value = 1 - sim$Py)

          set(sim, j = 'prodd0', value = 1 - sim$Pd)

          set(sim, j = 'poprisk', value = sim$prodp1)

        } else {

          # Calculate prodp1 as product of previous values
          set(sim, j = 'prodp1', value = sim$Py * (1 - sim$Pd) * sims[[i - 1]]$prodp0 * sims[[i - 1]]$prodd0)

          set(sim, j = 'poprisk', value = sims[[i - 1]]$poprisk + sim$prodp1)

          # Calculate probability of survival or D = 1
          set(sim, j = 'prodp0', value = (1 - sim$Py) * sims[[i - 1]]$prodp0)

          # Calculate probability of survival from D
          set(sim, j = 'prodd0', value = (1 - sim$Pd) * sims[[i - 1]]$prodd0)
        }
      }
      sims[[i]] <- sim
    }

    # efficiently bind sims
    sims <- rbindlist(sims)

    # initialize weights (only observed are weighted)
    sims[, .w := 1]

    # sort in proper order
    setorderv(sims, c(id, '.sid', time, '.sim'))

    if (prediction) {
      # calculate covariate means
      covs <- calc_obs_cov_means(covariate_fit, sims, id, time, by = c(id))

      # calculate predictions
      bylist <- c(id, time)
      if (outcome_fit$type == 'survival') {
        means <- sims[, .(poprisk = mean(poprisk)), by = bylist]
      } else if (outcome_fit$type == 'continuous') {
        means <- sims[, .(Ey = mean(Ey)), by = bylist]
      } else if (outcome_fit$type == 'binomial') {
        means <- sims[, .(Py = mean(Py)), by = bylist]
      }

      if (last_only) {
        rows <- which(means[[time]] == stop_time)
        means <- means[rows, ]
      }
    } else {

      # calculate covariate means
      covs <- calc_obs_cov_means(covariate_fit, sims, id, time)

      # calculate outcome means
      if (outcome_fit$type == 'survival') {
        means <- sims[, .(poprisk = mean(poprisk)), by = time]
      } else if (outcome_fit$type == 'continuous') {
        means <- sims[, .(Ey = mean(Ey)), by = time]
        rows <- which(means[[time]] == stop_time)
        means <- means[rows, ]
      } else if (outcome_fit$type == 'binomial') {
        means <- sims[, .(Py = mean(Py)), by = time]
        rows <- which(means[[time]] == stop_time)
        means <- means[rows, ]
      }

    }

    # flagged <- unique(
    #   paste(
    #     sims[is.na(sims[['treat']]) | is.na(sims[["cov1"]] | is.na(sims[["cov2"]])), ][['id_num']],
    #     sims[is.na(sims[['treat']]) | is.na(sims[["cov1"]] | is.na(sims[["cov2"]])), ][['.sim']],
    #     sep = "_"
    #   )
    # )
    # print(sims[paste(sims$id_num, sims$.sim, sep = "_") %in% flagged, ])
    if (!return_sims) {
      sims <- NULL
    }

    return(
      list(
        sims = sims,
        means = means,
        covs = covs
      )
    )

  }

#' Title
#'
#' @param uids
#' @param ind
#' @param gformula
#' @param splits
#' @param newdata
#' @param interventions
#' @param n_sims
#' @param n_samples
#' @param start_time
#' @param stop_time
#' @param bound_sims
#' @param last_only
#' @param prediction
#'
#' @return
#'
#' @examples
bootstrap_simulations <- function (uids,
                                   ind,
                                   gformula,
                                   splits,
                                   newdata,
                                   interventions,
                                   n_sims,
                                   n_samples,
                                   start_time,
                                   stop_time,
                                   bound_sims,
                                   last_only,
                                   prediction) {


  # bootstrap dataset
  blocks <- splits[ind]
  boot <- gformula$data[unlist(blocks), ]

  # re-fit models
  bgf <- gformula(
    outcome_model = gformula$outcome_model,
    covariate_model = gformula$covariate_model,
    compevent_model = gformula$compevent_model,
    censor_model = gformula$censor_model,
    visit_model = gformula$visit_model,
    data = boot,
    id = gformula$id,
    time = gformula$time,
    start_time = gformula$start_time,
    stop_time = gformula$stop_time,
    survival = gformula$survival,
    covs_baseline = gformula$covs_baseline,
    covs_tv = gformula$covs_tv
  )

  # re-simulate interventions
  sim_results <- run_gformula(
    gformula = bgf,
    newdata = if (!prediction) {
      if (is.null(newdata)) NULL else newdata[ind]
    } else {
      if (is.null(newdata)) gformula$data else newdata
    },
    interventions = interventions,
    n_sims = n_sims,
    n_samples = n_samples,
    n_boots = NULL,
    start_time = start_time,
    stop_time = stop_time,
    bound_sims = bound_sims,
    return_sims = FALSE,
    natural_course = FALSE,
    last_only = last_only
  )

  if (gformula$fit$outcome$type == "survival") {
    target <- "poprisk"
  } else if (gformula$fit$outcome$type == "continuous") {
    target <- "Ey"
  } else {
    target <- "Py"
  }

  statistic <- collapse_mean_covs(sim_results$means, sim_results$covs, target, prediction)

  return(statistic)
}
