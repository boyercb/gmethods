make_histories <- function(data, id, variables) {
  cumsum_vars <- gsub("^lag[0-9]+_", "", variables)
  cumsum_vars <- unique(cumsum_vars[grepl("cumsum_", cumsum_vars)])
  cumsum_targets <- gsub("cumsum_", "", cumsum_vars)

  cumavg_vars <- gsub("^lag[0-9]+_", "", variables)
  cumavg_vars <- unique(cumavg_vars[grepl("cumavg_", cumavg_vars)])
  cumavg_targets <- gsub("cumavg_", "", cumavg_vars)

  lagged_vars <- unique(variables[grepl("^lag[0-9]+_", variables)])
  lagged_targets <- gsub("^lag[0-9]+_", "", lagged_vars)
  lagged_vals <- gsub("^lag", "", unlist(regmatches(lagged_vars, regexec("^lag[0-9]+", lagged_vars))))
  lagged_vals <- as.list(as.numeric(lagged_vals))

  ts_vars <- gsub("^lag[0-9]+_", "", variables)
  ts_vars <- unique(ts_vars[grepl("ts_", ts_vars)])
  ts_targets <- gsub("ts_", "", ts_vars)

  # TODO: enforce having all lags below maximum encountered
  if (length(cumsum_vars) > 0) {
    data[, (cumsum_vars) := lapply(.SD, cumsum), by = id, .SDcols = cumsum_targets]
  }

  if (length(cumavg_vars) > 0) {
    data[, (cumavg_vars) := lapply(.SD, function(x)
        cumsum(x) / (1:.N)), by = id, .SDcols = cumavg_targets]
  }

  if (length(ts_vars) > 0) {
    data[, (ts_vars) := lapply(.SD, cumsum), by = id, .SDcols = ts_targets]

    for (j in seq_along(ts_vars)) {
      data[, (ts_vars[j]) := 1:.N - 1, by = c(id, ts_vars[j])]
    }
  }

  if (length(lagged_vars) > 0) {
    factors <- sapply(data[, .SD, .SDcols = lagged_targets], is.factor)
    if (any(factors)) {
      data[, (lagged_vars[!factors]) := mapply(data.table::shift,
                                               x = .SD,
                                               n = lagged_vals[!factors],
                                               fill = 0,
                                               SIMPLIFY = FALSE),
           by = id, .SDcols = lagged_targets[!factors]]

      reflevels <- lapply(lagged_targets[factors], function(x) levels(data[[x]])[1])

      for (j in seq_along(lagged_vars[factors])) {
        data[, (lagged_vars[factors][j]) := mapply(data.table::shift,
                                                   x = .SD,
                                                   n = lagged_vals[factors][j],
                                                   fill = reflevels[[j]],
                                                   SIMPLIFY = FALSE),
             by = id, .SDcols = lagged_targets[factors][j]]
      }

    } else {
      data[, (lagged_vars) := mapply(data.table::shift,
                                     x = .SD,
                                     n = lagged_vals,
                                     fill = 0,
                                     SIMPLIFY = FALSE),
          by = id, .SDcols = lagged_targets]
    }
  }

  return(data[])
}

update_cumulative_histories <- function(sim, sims, variable) {
  target <- gsub("(^cumsum_)|(^cumavg_)", "", variable)
  histories <- sapply(sims, function(s) s[[target]], simplify = TRUE)

  if (length(dim(histories)) >= 2L) {
    cumulative <- rowSums(histories, na.rm = TRUE)
  } else {
    cumulative <- histories
  }

  if (grepl("^cumsum_", variable)) {
    set(sim, j = variable, value = sim[[target]] + cumulative)
  } else {
    set(sim, j = variable, value = (sim[[target]] + cumulative) / (length(sims) + 1))
  }
}

update_ts_histories <- function(data, variable) {

  # extract lagged variables from histories
  target <- gsub("^ts_", "", variable)

  # update
  set(data,
      j = variable,
      value = ifelse(data[[target]] == 1, 0, data[[variable]] + 1)
  )

}

update_lagged_histories <- function(data, histories) {

  # extract lagged variables from histories
  vars <- histories[grepl("^lag[0-9]+_", histories)]

  # extract lagged values
  vals <- gsub("^lag", "", unlist(regmatches(vars, regexec("^lag[0-9]+", vars))))
  vals <- as.numeric(vals)

  # loop over variables and number of lags
  for (i in order(vals, decreasing = TRUE)) {
    if (vals[[i]] > 1) {
      set(data, j = vars[[i]],
          value = data[[gsub(vals[[i]], vals[[i]] - 1, vars[[i]])]])
    } else {
      set(data, j = vars[[i]],
          value = data[[gsub("^lag[0-9]+_", "", vars[[i]])]])
    }
  }

}


colSds <- function(x, na.rm = TRUE, tol = 1e-9) {
  if (na.rm) {
    n <- colSums(!is.na(x))
  } else {
    n <- nrow(x)
  }

  colVar <- colMeans(x * x, na.rm = na.rm) - (colMeans(x, na.rm = na.rm))^2
  colVar[abs(colVar) < tol] <- 0

  return(sqrt(colVar * n / (n - 1)))
}


get_header <- function(ints, ref) {
  header <- paste0("Parametric G-formula\n",
                   "Simulated risk under hypothetical interventions\n\n",
                   "Intervention \t Description\n")

  if (!is.null(ints)){
    for (i in seq_along(ints)){
      if ("Natural course" %in% ints) {
        j <- i - 1
      } else {
        j <- i
      }
      if (i == ref) {
        header <- paste0(header, j, "\t\t ", ints[i], " (*) \n")
      } else {
        header <- paste0(header, j, "\t\t ", ints[i], "\n")
      }
    }
  }

  header <- paste0(header, "\n* Reference intervention \n")

  header
}

get_footer <- function(size, n_samples, n_boots) {
  footer <- paste0("\nSample size = ", size,
                   "\nMonte Carlo samples = ", n_samples, "\n",
                   "Bootstrap samples = ", n_boots, "\n")

  footer
}

get_bootstrap_ci <- function(t, t0, conf_level, conf_type) {

  # calculate confidence intervals
  if (grepl("^norm", conf_type)) {
    if (length(dim(t)) < 2L) {
      se <- sd(t, na.rm = TRUE)
    } else {
      se <- colSds(t, na.rm = TRUE)
    }

    estimate <- t0
    conf.low <- estimate - se * qnorm((1 + conf_level) / 2)
    conf.high <- estimate + se * qnorm((1 + conf_level) / 2)

  } else if (grepl("^perc", conf_type)) {
    if (length(dim(t)) < 2L) {
      ll <- quantile(t, probs = (1 - conf_level) / 2, na.rm = TRUE)
      ul <- quantile(t, probs = (1 + conf_level) / 2, na.rm = TRUE)
    } else {
      ll <-
        apply(t, 2, function(x)
          quantile(x, probs = (1 - conf_level) / 2, na.rm = TRUE))
      ul <- apply(t, 2, function(x)
        quantile(x, probs = (1 + conf_level) / 2, na.rm = TRUE))
    }

    estimate <- t0
    conf.low <- ll
    conf.high <- ul

  } else if (grepl("^bias", conf_type)) {
    if (length(dim(t)) < 2L) {
      se <- sd(t, na.rm = TRUE)
      bias <- mean(t, na.rm = TRUE) - t0
    } else {
      se <- colSds(t, na.rm = TRUE)
      bias <- colMeans(t, na.rm = TRUE) - t0
    }

    estimate <- t0 - bias
    conf.low <- estimate - se * qnorm((1 + conf_level) / 2)
    conf.high <- estimate + se * qnorm((1 + conf_level) / 2)

  }

  return(cbind(estimate, conf.low, conf.high))
}

calc_np_means <-
  function(outcome_fit,
           covariate_fit,
           compevent_fit = NULL,
           censor_fit = NULL,
           visit_fit = NULL,
           data,
           id,
           time,
           censor_compevent = TRUE,
           stop_time = NULL) {

    dt <- data.table(data)

    # if there's censoring
    if (!is.null(censor_fit)) {
      censor <- as.character(formula(censor_fit)[2])

      # construct inverse probability of censoring weights
      dt[, .pC := 1 - predict(censor_fit, type = "response")]
      dt[, .pC := cumprod(.pC), by = c(id)]
      dt[, .w := 0][dt[[censor]] == 0, .w := 1 / .pC]

      # if there's a competing event that's also censored
      if (!is.null(compevent_fit) & outcome_fit$type == "survival") {
        compevent <- as.character(formula(compevent_fit)[2])

        if (censor_compevent) {
          # construct inverse probability of censoring weights for competing event
          dt[, .pD := 1 - predict(compevent_fit, type = "response")]
          dt[, .pD := cumprod(.pD), by = c(id)]
          dt[, .wD := 0][dt[[compevent]] == 0 | is.na(dt[[compevent]]), .wD := 1 / .pD]
          dt[, .w := .w + .wD]
        }
      }

      # if (trim_weights) {
      #
      # }

    } else {
      # set weights to one for everyone
      dt[, .w := 1]

    }

    obs_covs_means <- calc_obs_cov_means(covariate_fit, dt, id, time)
    obs_means <- calc_obs_means(outcome_fit, compevent_fit, dt, id, time, censor_compevent, stop_time)

    return(
      list(obs_covs_means = obs_covs_means, obs_means = obs_means)
    )
  }

calc_obs_cov_means <- function(covariate_fit, data, id, time, by = NULL) {
  covnames <- names(covariate_fit)
  types <- lapply(covariate_fit, get, x = "type")

  if (!is.null(by)) {
    bylist <- c(by, time)
  } else {
    bylist <- c(time)
  }

  if (any(types %in% "categorical")) {
    noncat <- which(!types %in% "categorical")
    cat <- which(types %in% "categorical")

    # for non-categorical types
    noncat_means <-
      data[, lapply(.SD, weighted.mean, w = .w), by = bylist, .SDcols = covnames[noncat]]

    # for categorical types
    # figure out how many levels
    lvls <- lapply(covnames[cat], function(v) levels(data[[v]]))
    lengths <- sapply(lvls, length)
    values <- unlist(lvls)

    SDnames <- rep(covnames[cat], each = lengths)
    newnames <- paste0(SDnames, "_", values)

    for (j in seq_along(SDnames)) {
      set(data, j = newnames[j], value = as.numeric(data[[SDnames[j]]] == values[j]))
    }

    cat_means <-
      data[, lapply(.SD, weighted.mean, w = .w), by = bylist, .SDcols = newnames]

    cov_means <- cat_means[noncat_means, on = time]
  } else {
    cov_means <-
      data[, lapply(.SD, weighted.mean, w = .w), by = bylist, .SDcols = covnames]
  }

  return(as.data.frame(cov_means))

}

calc_obs_means <- function(outcome_fit,
                           compevent_fit = NULL,
                           data,
                           id,
                           time,
                           censor_compevent,
                           stop_time = NULL) {

  outcome <- as.character(formula(outcome_fit)[2])

  if (outcome_fit$type == "survival") {

    if (is.null(compevent_fit) | censor_compevent) {
      out_means <- data[,
                        .(poprisk = sapply(.SD, weighted.mean, w = .w, na.rm = TRUE)),
                        by = time, .SDcols = outcome]
      out_means[, poprisk := poprisk * shift(cumprod(1 - poprisk), 1, 1)]
    } else {
      compevent <- as.character(formula(compevent_fit)[2])

      out_means <- data[,
                        lapply(.SD, weighted.mean, w = .w, na.rm = TRUE),
                        by = time, .SDcols = c(outcome, compevent)]

      setnames(out_means, c(time, 'poprisk', 'pD'))

      out_means[, poprisk := poprisk * (1 - pD) * shift(cumprod((1 - poprisk) * (1 - pD)), 1, 1)]
      out_means[, pD := NULL]
    }
    out_means[, poprisk := cumsum(poprisk)]

  } else if (outcome_fit$type == "continuous") {
    out_means <-
      data[, .(Ey = weighted.mean(.SD[[1]], .SD[[2]], na.rm = TRUE)),
           .SDcols = c(outcome, ".w")]
  } else if (outcome_fit$type == "binomial") {
    out_means <-
      data[, .(Py = weighted.mean(.SD[[1]], .SD[[2]], na.rm = TRUE)),
           .SDcols = c(outcome, ".w")]
  }

  return(as.data.frame(out_means))
}

# this function collapses simulated mean and covariates into a single vector
# for use during bootstrapping
collapse_mean_covs <- function(means, covs, target, prediction) {
  if (prediction) {
    remove_cols <- c(1,2)
  } else {
    remove_cols <- 1
  }
  l <- lapply(
    seq_along(means),
    function(i) {
      mean <- means[[i]][[target]]
      cov <- unlist(covs[[i]][, -remove_cols])
      return(c(mean, cov))
    })

  unlist(l)
}

# this function constructs a key to be used later for reconstructing mean and
# covs data.tables from collapsed vector
construct_key <- function(means, covs, target, interventions, prediction) {

  if (prediction) {
    remove_cols <- c(1,2)
  } else {
    remove_cols <- 1
  }

  last_idx <- 0

  # initialize key list
  key <- vector("list", length(means))

  for (i in seq_along(means)) {
    cov <- covs[[i]][, -remove_cols]

    # calculate length of means and covs
    n_means <- length(means[[i]][[target]])
    n_covs <- length(unlist(cov))

    # figure out corresponding indices in collapsed vector
    idx_means <- 1:n_means + last_idx

    # update last index placeholders
    last_cov_idx <- n_means + last_idx
    last_idx <- n_means + n_covs + last_idx

    # initialize covs indices list
    idx_covs <- vector("list", ncol(cov))

    # loop over each column and get corresponding indices in collapsed vector
    for (j in 1:ncol(cov)) {
      idx_covs[[j]] <- 1:length(cov[, j]) + last_cov_idx
      last_cov_idx <- length(cov[, j]) + last_cov_idx
    }

    # add names
    names(idx_covs) <- colnames(cov)

    # update key
    key[[i]]$means <- idx_means
    key[[i]]$covs <- idx_covs
  }

  # use intervention names
  names(key) <- names(interventions)

  return(key)
}

# use key to reconstruct mean data.tables from collapsed vector

use_means_key <- function(key, old, new, prediction) {
  key <- lapply(key, get, x = "means")

  means <- lapply(
    seq_along(key),
    function(i) {
      # get previous column names
      oldnames <- colnames(old[[i]])

      # extract data to be bound
      if (prediction) {
        df <- old[[i]][, 1:2]
      } else {
        df <- old[[i]][, 1]
      }

      # if coming from bootstrap with uncertainty
      if (length(dim(new)) > 1L) {
        # for single values (i.e. observed) combine rowwise
        if (length(dim(new[key[[i]], ])) > 1L) {
          df <- cbind(df, new[key[[i]], ])
        } else {
          df <- cbind(df, t(new[key[[i]], ]))
        }
        colnames(df) <- c(oldnames, "conf.low", "conf.high")
      } else {
        df <- cbind(df, new[key[[i]]])
        colnames(df) <- oldnames
      }
      return(df)
    })

  names(means) <- names(old)
  return(means)
}

use_covs_key <- function(key, old, new, prediction) {
  key <- lapply(key, get, x = "covs")

  covs <- lapply(
    seq_along(key),
    function(i) {

      # extract data to be bound
      if (prediction) {
        df <- old[[i]][, 1:2]

        # get previous column names
        oldnames <- colnames(old[[i]])[1:2]
      } else {
        df <- old[[i]][, 1]

        # get previous column names
        oldnames <- colnames(old[[i]])[1]
      }

      for (j in seq_along(key[[i]])) {
        covname <- names(key[[i]])[[j]]

        if (length(dim(new)) > 1L) {
          # for single values (i.e. observed) combine rowwise
          # coming from bootstrap with uncertainty
          df <- cbind(df, new[key[[i]][[j]], ])

          colnames(df) <-
            c(oldnames,
              covname,
              paste0(covname, ".low"),
              paste0(covname, ".high"))

        } else {
          df <- cbind(df, new[key[[i]][[j]]])

          colnames(df) <- c(oldnames, covname)
        }
        rownames(df) <- NULL
        oldnames <- colnames(df)

      }

      return(df)
    })

  names(covs) <- names(old)
  return(covs)
}

use_means_key_b <- function(key, old, new, prediction) {
  key <- lapply(key, get, x = "means")

  means <- lapply(
    seq_along(key),
    function(i) {
      df <- new[, key[[i]]]
      if (length(dim(df)) < 2L) {
        df <- matrix(df, ncol = 1)
      }
      if (prediction) {
        colnames(df) <- as.vector(paste(old[[i]][[1]], old[[i]][[2]], sep = "_"))
      } else {
        colnames(df) <- as.vector(old[[i]][[1]])
      }
      return(df)
    })

  names(means) <- names(old)
  return(means)
}

use_covs_key_b <- function(key, old, new, prediction) {
  key <- lapply(key, get, x = "covs")

  covs <- lapply(
    seq_along(key),
    function(i) {
      # loop over covs
      stats <- lapply(
        seq_along(key[[i]]),
        function(j) {
          # extract boot statistics
          df <- new[, key[[i]][[j]]]
          if (prediction) {
            colnames(df) <- as.vector(paste(old[[i]][[1]], old[[i]][[2]], sep = "_"))
          } else {
            colnames(df) <- as.vector(old[[i]][[1]])
          }
          return(df)
        })

      names(stats) <- names(key[[i]])

      return(stats)
    })

  names(covs) <- names(old)
  return(covs)
}
