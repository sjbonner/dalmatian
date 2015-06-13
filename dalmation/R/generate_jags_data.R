generateJAGSdata <- function(data,mean.model,variance.model,response=NULL,lower=NULL,upper=NULL){
    ## Generate data list
    jags.data <- list(n=nrow(data))
    
    ## Construct design matrices for mean model
    jags.data$mean.fixed <- model.matrix(mean.model$fixed$formula,data)
    jags.data[[paste0(mean.model$fixed$name,".n")]] <- ncol(jags.data$mean.fixed)

    if(!is.null(mean.model$random)){
        jags.data$mean.random <- model.matrix(mean.model$random$formula,data)
        jags.data[[paste0(mean.model$random$name,".n")]] <- ncol(jags.data$mean.random)
    }

    ## Construct design matrices for variance model
    jags.data$variance.fixed <- model.matrix(variance.model$fixed$formula,data)
    jags.data[[paste0(variance.model$fixed$name,".n")]] <- ncol(jags.data$variance.fixed)

    if(!is.null(variance.model$random)){
        jags.data$variance.random <- model.matrix(variance.model$random$formula,data)
        jags.data[[paste0(variance.model$random$name,".n")]] <- ncol(jags.data$variance.random)
    }

    ## Add response variable
    if(!is.null(response))
        jags.data$y <- data[,response]
    else if(!is.null(lower) && !is.null(upper)){
        jags.data$lower <- data[,lower]
        jags.data$upper <- data[,upper]
        jags.data$dummy <- rep(1,nrow(data))
    }
    else
        stop("You must either specify the exact response value or both lower and upper bounds for rounding.\n\n")

    ## Return data list
    jags.data
}

