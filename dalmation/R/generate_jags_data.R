generateJAGSdata <- function(df,mean.model,variance.model,response=NULL,lower=NULL,upper=NULL){
    ## Generate data list
    jags.data <- list(n=nrow(df))
    
    ## Construct design matrices for mean model
    jags.data$mean.fixed <- model.matrix(mean.model$fixed$formula,df)
    jags.data[[paste0(mean.model$fixed$name,".n")]] <- ncol(jags.data$mean.fixed)

    if(!is.null(mean.model$random)){
        jags.data$mean.random <- model.matrix(mean.model$random$formula,df)

        tmp <- attr(jags.data$mean.random,"assign")

        jags.data[[paste0(mean.model$random$name,".ncomponents")]] <- max(tmp)
        jags.data[[paste0(mean.model$random$name,".neffects")]] <- length(tmp)
        jags.data[[paste0(mean.model$random$name,".levels")]] <- tmp
    }

    ## Construct design matrices for variance model
    jags.data$variance.fixed <- model.matrix(variance.model$fixed$formula,df)
    jags.data[[paste0(variance.model$fixed$name,".n")]] <- ncol(jags.data$variance.fixed)

    if(!is.null(variance.model$random)){
        jags.data$variance.random <- model.matrix(variance.model$random$formula,df)

        tmp <- attr(jags.data$variance.random,"assign")

        jags.data[[paste0(variance.model$random$name,".ncomponents")]] <- max(tmp)
        jags.data[[paste0(variance.model$random$name,".neffects")]] <- length(tmp)
        jags.data[[paste0(variance.model$random$name,".levels")]] <- tmp
    }

    ## Add response variable
    if(!is.null(response))
        jags.data$y <- df[,response]
    else if(!is.null(lower) && !is.null(upper)){
        jags.data$lower <- df[,lower]
        jags.data$upper <- df[,upper]
        jags.data$dummy <- rep(1,nrow(df))
    }
    else
        stop("You must either specify the exact response value or both lower and upper bounds for rounding.\n\n")

    ## Return data list
    jags.data
}

