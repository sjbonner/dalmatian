generateJAGSdata <-
  function(df,
           family,
           mean.model,
           dispersion.model,
           joint.model,
           response = NULL,
           ntrials = NULL,
           lower = NULL,
           upper = NULL,
           include.checks = TRUE,
           drop.levels = TRUE,
           drop.missing = TRUE) {
    ## Remove unused factor level from data frame
    if (drop.levels){
      df <- droplevels(df)
    }

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

    ## Construct design matrices for dispersion model
    jags.data$dispersion.fixed <-
      model.matrix(dispersion.model$fixed$formula, df)
    jags.data[[paste0(dispersion.model$fixed$name, ".n")]] <-
      ncol(jags.data$dispersion.fixed)

    if (!is.null(dispersion.model$random)) {
      jags.data$dispersion.random <-
        model.matrix(dispersion.model$random$formula, df)

      tmp <- attr(jags.data$dispersion.random, "assign")

      jags.data[[paste0(dispersion.model$random$name, ".ncomponents")]] <-
        max(tmp)
      jags.data[[paste0(dispersion.model$random$name, ".neffects")]] <-
        length(tmp)
      jags.data[[paste0(dispersion.model$random$name, ".levels")]] <-
        tmp
    }
    
    ## Construct design matrices for joint components of the model
    if(!is.null(joint.model)){
      if(!is.null(joint.model$fixed)){
        jags.data$joint.fixed <-
          model.matrix(joint.model$fixed$formula, df)
        jags.data[[paste0(joint.model$fixed$name, ".n")]] <-
          ncol(jags.data$joint.fixed)
      }

      if (!is.null(joint.model$random)) {
        jags.data$joint.random <-
          model.matrix(joint.model$random$formula, df)
        
        tmp <- attr(jags.data$joint.random, "assign")
        
        jags.data[[paste0(joint.model$random$name, ".ncomponents")]] <-
          max(tmp)
        jags.data[[paste0(joint.model$random$name, ".neffects")]] <-
          length(tmp)
        jags.data[[paste0(joint.model$random$name, ".levels")]] <-
          tmp
      }
    }

    ## Construct weights for dispersion model
    if (!is.null(dispersion.model$weights))
      jags.data$weights <- df[[dispersion.model$weights]]

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

    ## Add number of trials for betabinomial
    if(family == "betabin"){
      if(is.null(ntrials))
        stop("You must specify the number of independent trials for each observation of the beta-binomial model.\n\n")
      
      jags.data$m <- dplyr::pull(df,ntrials)
    }

    if(include.checks){
      ## Add dummy variables for checking support of mean and dispersion
      jags.data$mean.check <- rep(1, nrow(df))
      jags.data$disp.check <- rep(1, nrow(df))
    }
    
    ## Return data list
    jags.data
  }
