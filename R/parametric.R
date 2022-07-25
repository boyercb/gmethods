#' Fit GLM on Covariate
#'
#' This internal function fits a generalized linear model (GLM) for a single covariate using the observed data.
#'
#' @param formula
#' @param family
#' @param link
#' @param data
#' @param control
#'
#' @return            Fitted model for the covariate.
#' @keywords internal

fit_glm <- function(formula, family, link = NULL, data, control = NULL) {
  if (is.null(link)) {
    famtext <- paste(family, "()", sep = "")
  } else {
    famtext <- paste(family, "(link = ", link, ")", sep = "")
  }

  # Fit GLM for covariate using user-specified formula
  if (!is.null(control)){
    fit <- stats::glm(
      formula,
      family = eval(parse(text = famtext)),
      data = data,
      y = TRUE,
      control = control
    )
  } else {
    fit <- stats::glm(
      formula,
      family = eval(parse(text = famtext)),
      data = data,
      y = TRUE
    )
  }

  fit$rmse <- add_rmse(fit)
  fit$stderrs <- add_stderr(fit)
  fit$vcov <- add_vcov(fit)
  fit$outcome <- all.vars(update(formula, . ~ 1))

  fit$type <- family

  if (family == "gaussian") {
    fit$type <- "normal"
  }

  return(fit)
}

#' Fit Multinomial Model on Covariate
#'
#' This internal function fits a multinomial regression model for a categorical covariate using the observed data.
#'
#' @param formula
#' @param data
#' @param control
#'
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal

fit_multinomial <- function(formula, data, control = NULL){
  if (!is.null(control)) {
    args <- c(list(formula = formula, data = data), trace = FALSE, control)
    fit <- do.call(nnet::multinom, args = args)
  } else {

    fit <- nnet::multinom(formula = formula,
                          data = data,
                          trace = FALSE)
  }

  fit$rmse <- add_rmse(fit)
  # fit$stderr <- add_stderr(fit)
  # fit$vcov <- add_vcov(fit)
  fit$outcome <- all.vars(update(formula, . ~ 1))
  fit$type <- 'categorical'

  return(fit)
}

#' Fit Zero-Inflated Normal Model on Covariate
#'
#' This internal function models a zero-inflated normal distribution through the combined
#' use of a generalized linear model (GLM) fit on a zero vs. non-zero indicator
#' and a GLM fit on all non-zero values.
#'
#' @param formula
#' @param family
#' @param link
#' @param data
#' @param control
#'
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal
#' @import data.table

fit_zeroinfl_normal <- function(formula, link = NULL, data, control = NULL) {
  if (is.na(link)) {
    famtext <- "gaussian()"
  } else {
    famtext <- paste("gaussian(link = ", link, ")", sep = "")
  }

  cov.name <- all.vars(update(formula, . ~ 1))

  # Take log to ensure that no negative values are predicted
  data[, paste("log_", cov.name, sep = "")] <- 0
  data[data[[cov.name]] != 0, ][[paste("log_", cov.name, sep = "")]] <-
    log(data[data[[cov.name]] != 0, ][[cov.name]])

  # Fit binomial model on indicator of whether covariate is 0 or not
  if (!is.null(control)) {
    fit0 <-
      stats::glm(stats::update(formula, stats::as.formula(paste0("I(", cov.name, " != 0) ~ ."))),
                 family = stats::binomial(),
                 control = control[[1]], data = data)
  } else {
    fit0 <-
      stats::glm(
        stats::update(formula, stats::as.formula(paste0("I(", cov.name, " != 0) ~ ."))),
        family = stats::binomial(),
        data = data
      )
  }

  # Fit Gaussian model on data points for which covariate does not equal 0
  if (!is.null(control)) {
    fit <- stats::glm(stats::update(formula, stats::as.formula(paste0("log_", cov.name, " ~ ."))),
                      family = eval(parse(text = famtext)),
                      control = control[[2]],
                      data = data[data[[cov.name]] != 0, ])
  } else {
    fit <- stats::glm(stats::update(formula, stats::as.formula(paste0("log_", cov.name, " ~ ."))),
                      family = eval(parse(text = famtext)),
                      data = data[data[[cov.name]] != 0, ])
  }

  fit0$rmse <- add_rmse(fit0)
  fit$rmse <- add_rmse(fit)

  fit0$stderr <- add_stderr(fit0)
  fit$stderr <- add_stderr(fit)

  fit0$vcov <- add_vcov(fit0)
  fit$vcov <- add_vcov(fit)

  fit0$outcome <- all.vars(update(formula, . ~ 1))
  fit$outcome <- all.vars(update(formula, . ~ 1))

  fit0$type <- 'zero-inflated normal'
  fit$type <- 'zero-inflated normal'

  return(list(
    fits = list(zero = fit0, normal = fit),
    type = 'zero-inflated normal',
    rmse = fit$rmse,
    outcome = fit$outcome
  ))
}

#' Fit Hurdle Model for Counts on Covariate
#'
#' This internal function models a covariate using a hurdle model, i.e. a
#' two-component model in which zeros are modeled seperately from a truncated
#' distribution for positive counts. The count model can either be a truncated
#' Poisson or negative binomial model. The zero or "hurdle" component is
#' modeled with a binomial GLM (e.g. logistic regression).
#'
#' @param formula
#' @param family
#' @param link
#' @param data
#' @param control
#'
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal
#' @import data.table
fit_hurdle <- function(formula, family, link, data, control = NULL) {
  # check that pscl is installed
  if (!requireNamespace("pscl", quietly = TRUE)) {
    stop("Package \"pscl\" must be installed to use this function.")
  }

  stopifnot(link %in% c("logit", "probit", "cloglog", "cauchit", "log"),
            "Link must be one of: \"logit\", \"probit\", \"cloglog\", \"cauchit\", \"log\"")

  # Fit GLM for covariate using user-specified formula
  if (!is.null(control)){
    fit <- pscl::hurdle(
      formula,
      data = data,
      dist = family,
      link = link,
      y = TRUE,
      control = control
    )
  } else {
    fit <- pscl::hurdle(
      formula,
      data = data,
      dist = family,
      link = link,
      y = TRUE
    )
  }

  fit$rmse <- add_rmse(fit)
  fit$stderrs <- add_stderr(fit)
  fit$vcov <- add_vcov(fit)
  fit$outcome <- all.vars(update(formula, . ~ 1))

  fit$type <- family

  return(fit)

}

#' Fit Zero-inflated Model for Counts on Covariate
#'
#' This internal function models a covariate using a zero-inflated count model,
#' i.e. a two-component mixture model in which observations are assumed to
#' arise from a mixture distribution composed of a 0/1 distribution such as the
#' bernoulli/binomial and a proper count distribution. It differs from a "hurdle"
#' model in that there are two sources of zeros: zeros may come from the 0/1
#' distribution or from the proper count distribution.
#'
#' @param formula
#' @param family
#' @param link
#' @param data
#' @param control
#'
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal
#' @import data.table
fit_zeroinfl_count <- function(formula, family, link, data, control = NULL) {
  # check that pscl is installed
  if (!requireNamespace("pscl", quietly = TRUE)) {
    stop("Package \"pscl\" must be installed to use this function.")
  }

  stopifnot(link %in% c("logit", "probit", "cloglog", "cauchit", "log"),
            "Link must be one of: \"logit\", \"probit\", \"cloglog\", \"cauchit\", \"log\"")
  # Fit GLM for covariate using user-specified formula
  if (!is.null(control)){
    fit <- pscl::zeroinfl(
      formula,
      data = data,
      dist = family,
      link = link,
      y = TRUE,
      control = control
    )
  } else {
    fit <- pscl::zeroinfl(
      formula,
      data = data,
      link = link,
      dist = family,
      y = TRUE
    )
  }

  fit$rmse <- add_rmse(fit)
  fit$stderrs <- add_stderr(fit)
  fit$vcov <- add_vcov(fit)
  fit$outcome <- all.vars(update(formula, . ~ 1))

  fit$type <- family

  return(fit)

}

#' Fit Bounded Normal Model on Covariate
#'
#' This internal function models a covariate using a "bounded normal" distribution
#' by first standardizing the covariate values to the range [0, 1], noninclusive,
#' then fitting a generalized linear model (GLM) under the Gaussian family function.
#'
#' @param formula
#' @param family
#' @param link
#' @param data
#' @param control
#'
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal
#' @import data.table

fit_bounded_continuous <- function(formula, link = NULL, data, control = NULL) {

  cov.name <- all.vars(update(formula, . ~ 1))

  data[, paste("norm_", cov.name, sep = "")] <-
    (data[[cov.name]] - min(data[[cov.name]])) /
    (max(data[[cov.name]]) - min(data[[cov.name]]))

  if (!is.null(link)) {
    if (!is.null(control)) {
      fit <- stats::glm(
        stats::update(formula, stats::as.formula(paste0("norm_", cov.name, " ~ ."))),
        family = stats::gaussian(link = link),
        data = data,
        y = TRUE,
        control = control
      )
    } else {
      fit <- stats::glm(
        stats::update(formula, stats::as.formula(paste0("norm_", cov.name, " ~ ."))),
        family = stats::gaussian(link = link),
        data = data,
        y = TRUE
      )
    }
  } else {
    if (!!is.null(control)) {
      fit <- stats::glm(
        stats::update(formula, stats::as.formula(paste0("norm_", cov.name, " ~ ."))),
        family = stats::gaussian(),
        data = data,
        y = TRUE,
        control = control
      )
    } else {
      fit <- stats::glm(
        stats::update(formula, stats::as.formula(paste0("norm_", cov.name, " ~ ."))),
        family = stats::gaussian(),
        data = data,
        y = TRUE
      )
    }
  }

  fit$rmse <- add_rmse(fit)
  fit$stderr <- add_stderr(fit)
  fit$vcov <- add_vcov(fit)
  fit$outcome <- cov.name
  fit$type <- 'bounded normal'
  return(fit)
}

#' Fit Truncated Normal Model on Covariate
#'
#' This internal function models a covariate using a normal distribution truncated on one
#' side at a user-specified cutoff.
#'
#' @param formula
#' @param family
#' @param link
#' @param data
#' @param control
#' @param direction
#' @param point
#'
#' @return            Fitted model for the covariate at index \eqn{j}.
#' @keywords internal

fit_trunc_normal <- function(formula, link = NULL, data, control = NULL, direction = NULL, point = NULL) {
  # check that truncreg is installed
  if (!requireNamespace("truncreg", quietly = TRUE)) {
    stop("Package \"truncreg\" must be installed to use this function.")
  }

  if (!is.null(control)) {
    args <- c(list(
      formula = formula,
      data = data,
      point = point,
      direction = direction,
      y = TRUE
    ),
    control)

    fit <- do.call(truncreg::truncreg, args = args)
  } else {
    fit <- truncreg::truncreg(
      formula,
      data = data,
      point = point,
      direction = direction,
      y = TRUE
    )
  }

  fit$rmse <- add_rmse(fit)
  fit$stderr <- add_stderr(fit)
  fit$vcov <- add_vcov(fit)
  fit$outcome <- all.vars(update(formula, . ~ 1))
  fit$type <- "truncated normal"
  return(fit)
}

#' Fit Covariate Models
#'
#' This internal function fits a model for each covariate using the observed data.
#'
#' @param models
#' @param data
#'
#' @return                A list of fitted models, one for each covariate in \code{covnames}.
#' @keywords internal
#' @import data.table

fit_covariate_model <- function(model, data, time = "time", start_time = 0, gam = FALSE) {

  # fit to all observations after baseline/start time
  subdata <- data[data[[time]] > start_time, ]

  if (!is.null(model$subset)) {
    subdata <- subset(subdata, eval(parse(text = model$subset)))
  }
  if (!is.null(model$restrict)) {
    subdata <- subset(subdata, eval(parse(text = model$restrict$subset)))
  }

  if (gam) {
    fit <- switch(
      model$family,
      'binomial' = gam::gam(
        formula = model$formula,
        family = binomial(link = model$link),
        data = subdata
      ),

      'normal' = gam::gam(
        formula = model$formula,
        family = gaussian,
        data = subdata
      )
    )
    fit$rmse <- add_rmse(fit)
    fit$stderrs <- add_stderr(fit)
    fit$vcov <- add_vcov(fit)
    fit$outcome <- all.vars(update(model$formula, . ~ 1))
    if (model$family == 'normal') {
      fit$type <- 'normal'
    } else if (model$family == 'binomial') {
      fit$type <- 'binomial'
    }

  } else {
    fit <- switch(
      model$family,
      'binomial' = fit_glm(
        formula = model$formula,
        family = 'binomial',
        link = model$link,
        data = subdata,
        control = model$control
      ),

      'normal' = fit_glm(
        formula = model$formula,
        family = 'gaussian',
        link = model$link,
        data = subdata,
        control = model$control
      ),

      'poisson' = fit_glm(
        formula = model$formula,
        family = 'poisson',
        link = model$link,
        data = subdata,
        control = model$control
      ),

      'categorical' = fit_multinomial(
        formula = model$formula,
        data = subdata,
        control = model$control
      ),

      'zero-inflated normal' = fit_zeroinfl_normal(
        formula = model$formula,
        link = model$link,
        data = subdata,
        control = model$control
      ),

      'poisson hurdle' = fit_hurdle_count(
        formula = model$formula,
        link = model$link,
        family = "poisson",
        data = subdata,
        control = model$control
      ),

      'negbin hurdle' = fit_hurdle_count(
        formula = model$formula,
        link = model$link,
        family = "negbin",
        data = subdata,
        control = model$control
      ),

      'zero-inflated poisson' = fit_zeroinfl_count(
        formula = model$formula,
        link = model$link,
        family = "poisson",
        data = subdata,
        control = model$control
      ),

      'zero-inflated negbin' = fit_zeroinfl_count(
        formula = model$formula,
        link = model$link,
        family = "negbin",
        data = subdata,
        control = model$control
      ),

      'bounded normal' = fit_bounded_continuous(
        formula = model$formula,
        link = model$link,
        data = subdata,
        control = model$control
      ),

      'truncated normal' = fit_trunc_normal(
        formula = model$formula,
        link = model$link,
        data = subdata,
        control = model$control,
        point = model$point,
        direction = model$direction
      )
    )
  }

  fit$restrict <- model$restrict

  return(fit)
}

#' Fit Outcome Model
#'
#' This internal function fits a generalized linear model (GLM) for the outcome variable using the observed data.
#'
#' @param data           A data frame, list or environment (or object coercible by
#'                       as.data.frame to a data frame) containing the variables in the model.
#'                       If not found in data, the variables are taken from environment(formula),
#'                       typically the environment from which glm is called.
#' @param EOF            Only use end of follow up time point when fitting outcome model
#' @param model          Model statement for the outcome variable.
#' @param time           Time variable.
#'
#' @return               Fitted model for the outcome variable.
#' @keywords internal
#' @import data.table

fit_outcome_model <- function(model, data, time = NULL, stop_time = NULL, survival = FALSE, gam = FALSE) {
  family <- switch(
    model$family,
    "normal" = stats::gaussian(link = model$link),
    "binomial" = stats::binomial(link = model$link)
  )

  outcome <- all.vars(update(model$formula, . ~ 1))

  if (!survival) {
    if (is.null(stop_time)) {
      stop("Must specify stop time for EOF.")
    }
    data <- data[data[[time]] == stop_time, ]
  }

  if (gam) {
    fit_function <- gam::gam
  } else {
    fit_function <- stats::glm
  }

  if (!is.null(model$subset)) {

    if (model$family == "normal") {
      data[, paste("norm_", outcome, sep = "")] <-
        (data[[outcome]] - min(data[[outcome]]))/(max(data[[outcome]]) - min(data[[outcome]]))

      fitY <- fit_function(
        stats::as.formula(paste("norm_", model, sep = "")),
        family = family,
        data = subset(data, eval(parse(text = model$subset)))
      )
    } else {
      fitY <- fit_function(
        formula = model$formula,
        family = family,
        data = subset(data, eval(parse(text = model$subset)))
      )
    }
  } else { # Case where there are no restrictions on outcome variable
    # Fit GLM for outcome variable using user-specified model and entire dataset
    fitY <- fit_function(
      formula = model$formula,
      family = family,
      data = data
    )
  }

  fitY$rmse <- add_rmse(fitY)
  fitY$stderr <- add_stderr(fitY)
  fitY$vcov <- add_vcov(fitY)
  fitY$outcome <- outcome

  if (!survival) {
    if (model$family == "normal") {
      fitY$type <- 'continuous'
    } else {
      fitY$type <- 'binomial'
    }
  } else {
    fitY$type <- 'survival'
  }

  return(fitY)
}

#' Fit Competing Event Model
#'
#' This internal function fits a generalized linear model (GLM) for the competing event variable using the observed data.
#'
#' @param data                   A data frame, list or environment (or object coercible by
#'                               as.data.frame to a data frame) containing the variables in the model.
#'                               If not found in data, the variables are taken from environment(formula),
#'                               typically the environment from which glm is called.
#' @param model                  Model statement for competing event variable.
#'
#' @return                       Fitted model for the competing event variable.
#' @keywords internal

fit_compevent_model <- function(model, data, use_gam = FALSE) {
  if (!model$family == "binomial") {
    stop("Currently only binomial models for competing events are allowed")
  }

  # if (use_gam) {
  #   fit_function <- gam::gam
  # } else {
  #   fit_function <- stats::glm
  # }
  fit_function <- stats::glm

  if (!is.null(model$subset)) {
    # Check for restrictions on compevent event variable modeling
    # Fit GLM assuming binomial distribution
    # Restrict data used for modeling to rows where condition is true
    fitD <- fit_function(
      formula = model$formula,
      family = stats::binomial(link = model$link),
      data = subset(data, eval(parse(text = model$subset)))
    )
  } else {
    # Case where there are no restrictions on competing event variable
    #Fit GLM assuming binomial distribution
    fitD <- fit_function(
      formula = model$formula,
      family = stats::binomial(link = model$link),
      data = data
    )
  }

  fitD$rmse <- add_rmse(fitD)
  fitD$stderr <- add_stderr(fitD)
  fitD$vcov <- add_vcov(fitD)
  fitD$outcome <- all.vars(update(model$formula, . ~ 1))
  fitD$type <- 'survival'
  return(fitD)
}

trim_glm <- function(fit) {
  fit$y <- c()
  fit$model <- c()

  fit$residuals <- c()
  fit$fitted.values <- c()
  fit$effects <- c()
  fit$qr$qr <- c()
  fit$linear.predictors <- c()
  fit$weights <- c()
  fit$prior.weights <- c()
  fit$data <- c()

  fit$family$variance <- c()
  fit$family$dev.resids <- c()
  fit$family$aic <- c()
  fit$family$validmu <- c()
  fit$family$simulate <- c()
  attr(fit$terms, ".Environment") <- c()
  attr(fit$formula, ".Environment") <- c()

  return(fit)
}

trim_truncreg <- function(fit) {
  fit$gradientObs <- c()
  fit$fitted.values <- c()
  fit$y <- c()

  return(fit)
}

trim_multinom <- function(fit) {
  fit$fitted.values <- c()
  fit$residuals <- c()
  fit$weights <- c()
  fit$y <- c()

  return(fit)
}

add_rmse <- function(fit) {
  return(sqrt(mean((fit$y - stats::fitted(fit))^2)))
}

add_grf_rmse <- function(rf) {
  if (rf$type %in% c("binomial grf", "survival")) {
    return(sqrt(mean((as.numeric(rf$Y.orig) - 1 - rf$predictions[, 2])^2)))
  } else {
    return(sqrt(mean((rf$Y.orig - rf$predictions)^2)))
  }
}

add_stderr <- function(fit) {
  if (any(class(fit) == 'multinom')) {
    return(summary(fit)$coefficients)
  } else {
    return (stats::coefficients(summary(fit))[, "Std. Error"])
  }
}

add_vcov <- function(fit) {
  return(stats::vcov(fit))
}


winsorize <- function(x, pct) {
  upr <- quantile(x, pct + (1 - pct) / 2, na.rm = TRUE)
  lwr <- quantile(x, (1 - pct) / 2, na.rm = TRUE)

  x[x > upr] <- upr
  x[x < lwr] <- lwr

  return(x)
}

draw_grf_binomial <- function(N, size, p.fit, data) {
  if (is.factor(p.fit$Y.orig)) {
    stats::rbinom(
      n = N,
      size = size,
      prob = stats::predict(p.fit, type = 'response', newdata = data)$predictions[,2]
    )
  } else {
    stats::rbinom(
      n = N,
      size = size,
      prob = stats::predict(p.fit, type = 'response', newdata = data)$predictions
    )
  }
}

draw_grf_normal <- function(N, mu.fit, sd.hat, data) {
  stats::rnorm(
    n = N,
    mean = stats::predict(mu.fit, type = 'response', newdata = data)$predictions,
    sd = sd.hat
  )
}

draw_binomial <- function(N, size, p.fit, data) {
  stats::rbinom(
    n = N,
    size = size,
    prob = stats::predict(p.fit, type = 'response', newdata = data)
  )
}

draw_normal <- function(N, mu.fit, sd.hat, data) {
  stats::rnorm(
    n = N,
    mean = stats::predict(mu.fit, type = 'response', newdata = data),
    sd = sd.hat
  )
}

draw_multinomial <- function(p.fit, data) {
  stats::predict(p.fit, type = 'class', newdata = data)
}

draw_poisson <- function(N, lambda.fit, data) {
  stats::rpois(
    n = N,
    lambda = predict(lambda.fit, type = 'response', newdata = data)
  )
}

draw_trunc_poisson <- function(N, lambda) {
  u <- runif(N, min = ppois(0, lambda), max = 1)
  qpois(u, lambda = lambda)
}

draw_trunc_negbin <- function(N, lambda, size) {
  u <- runif(N, min = pnbinom(0, mu = lambda, size = size), max = 1)
  qnbinom(u, mu = lambda, size = size)
}

draw_zeroinfl_poisson <- function(N, size, fit, data) {
  nonzero <- draw_binomial(
    N = N,
    size = size,
    p.fit = pscl:::predict.hurdle(fit, type = "zero", newdata = data),
    data = data
  )

  response <- draw_poisson(
    N = N,
    lambda.fit = pscl:::predict.hurdle(fit, type = "count"),
    data = data
  )

  nonzero * response
}

draw_zeroinfl_negbin <- function(N, size, fit, data) {
  nonzero <- draw_binomial(
    N = N,
    size = size,
    p.fit = pscl:::predict.hurdle(fit, type = "zero", newdata = data),
    data = data
  )

  response <- rnbinom(
    n = N,
    mu = pscl:::predict.hurdle(fit, type = "count", newdata = data),
    size = fit$theta
  )

  nonzero * response
}

draw_poisson_hurdle <- function(N, size, fit, data) {
  nonzero <- draw_binomial(
    N = N,
    size = size,
    p.fit = pscl:::predict.hurdle(fit, type = "zero", newdata = data),
    data = data
  )

  response <- draw_trunc_poisson(
    N = N,
    lambda = pscl:::predict.hurdle(fit, type = "count", newdata = data)
  )

  nonzero * response
}

draw_negbin_hurdle <- function(N, size, fit, data) {
  nonzero <- draw_binomial(
    N = N,
    size = size,
    p.fit = pscl:::predict.hurdle(fit, type = "zero", newdata = data),
    data = data
  )

  response <- draw_trunc_negbin(
    N = N,
    lambda = pscl:::predict.hurdle(fit, type = "count", newdata = data),
    size = fit$theta
  )

  nonzero * response
}

draw_zero_inflated_normal <- function(N, size, p.fit, mu.fit, sd.hat, data) {
  nonzero <- draw_binomial(
    N = N,
    size = size,
    p.fit = p.fit,
    data = data
  )

  response <- draw_normal(
    N = N,
    mu.fit = mu.fit,
    sd.hat = sd.hat,
    data = data
  )

  nonzero * exp(response)
}

draw_bounded_normal <- function(N, mu.fit, sd.hat, range, data) {
  response <- draw_normal(
    N = N,
    mu.fit = mu.fit,
    sd.hat = sd.hat,
    data = data
  )

  (response * (range[2] - range[1])) + range[1]
}

draw_trunc_normal <- function(N, mu.fit, sd.hat, direction, point, data) {
  if (direction == 'left') {
    a <- point
    b <- Inf
  } else if (direction == 'right') {
    a <- -Inf
    b <- point
  }

  truncnorm::rtruncnorm(
    n = N,
    mean = predict(mu.fit, type = 'response', newdata = data),
    sd = sd.hat,
    a = a,
    b = b
  )
}

draw_custom <- function(N, mu.fit, sd.hat, data, pred.fun = stats::predict()) {
  stats::rnorm(
    n = N,
    mean = stats::predict(mu.fit, type = 'response', newdata = data),
    sd = rmse.hat
  )
}


