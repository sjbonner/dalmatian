##' @importFrom stats model.matrix
generateJAGSdata <-
  function(df,
           mean.model,
           variance.model,
           response = NULL,
           lower = NULL,
           upper = NULL,
           drop.levels = TRUE,
           drop.missing = TRUE) {
    ## Remove unused factor level from data frame
    if (drop.levels)
      df <- gdata::drop.levels(df, reorder = FALSE)

    ## Remove records with missing responses
    if (is.null(response)) {
      miss <- intersect(which(!is.finite(df[[upper]])),
                        which(!is.finite(df[[lower]])))
    }
    else if (drop.missing)
      miss <- which(!is.finite(df[[response]]))

    if (length(miss) > 0) {
      df <- df[-miss,]

      warning("Removed ",
              length(miss),
              " observations with missing response values.\n")
    }

    ## Initialize data list
    jags.data <- list(n = nrow(df))

    ## Construct design matrices for mean model
    jags.data$mean.fixed <-
      model.matrix(mean.model$fixed$formula, df)
    jags.data[[paste0(mean.model$fixed$name, ".n")]] <-
      ncol(jags.data$mean.fixed)

    if (!is.null(mean.model$random)) {
      jags.data$mean.random <- model.matrix(mean.model$random$formula, df)

      tmp <- attr(jags.data$mean.random, "assign")

      jags.data[[paste0(mean.model$random$name, ".ncomponents")]] <-
        max(tmp)
      jags.data[[paste0(mean.model$random$name, ".neffects")]] <-
        length(tmp)
      jags.data[[paste0(mean.model$random$name, ".levels")]] <-
        tmp

      if (!is.null(mean.model$random$sigma))
        jags.data[[paste0("sd.", mean.model$random$name)]] <-
        mean.model$random$sigma ^ 2
    }

    ## Construct design matrices for variance model
    jags.data$variance.fixed <-
      model.matrix(variance.model$fixed$formula, df)
    jags.data[[paste0(variance.model$fixed$name, ".n")]] <-
      ncol(jags.data$variance.fixed)

    if (!is.null(variance.model$random)) {
      jags.data$variance.random <-
        model.matrix(variance.model$random$formula, df)

      tmp <- attr(jags.data$variance.random, "assign")

      jags.data[[paste0(variance.model$random$name, ".ncomponents")]] <-
        max(tmp)
      jags.data[[paste0(variance.model$random$name, ".neffects")]] <-
        length(tmp)
      jags.data[[paste0(variance.model$random$name, ".levels")]] <-
        tmp
    }

    ## Construct weights for variance model
    if (!is.null(variance.model$weights))
      jags.data$weights <- df[[variance.model$weights]]

    ## Add response variable
    if(("tbl" %in% class(df))){
      if (!is.null(response))
        jags.data$y <- dplyr::pull(df,response)
      else if (!is.null(lower) && !is.null(upper)) {
        jags.data$lower <- dplyr::pull(df,lower)
        jags.data$upper <- dplyr::pull(df,upper)
        jags.data$dummy <- rep(1, nrow(df))
      }
      else
        stop(
          "You must either specify the exact response value or both lower and upper bounds for rounding.\n\n"
        )
    }
    else{
      if (!is.null(response))
        jags.data$y <- df[, response]
      else if (!is.null(lower) && !is.null(upper)) {
        jags.data$lower <- df[, lower]
        jags.data$upper <- df[, upper]
        jags.data$dummy <- rep(1, nrow(df))
      }
      else
        stop(
          "You must either specify the exact response value or both lower and upper bounds for rounding.\n\n"
        )
    }

    ## Return data list
    jags.data
  }
