generateJAGSinits <- function(mean.model,variance.model,jags.data){

    ## Initialize list
    inits <- list()
    
    ## Generate initial values for parameters of mean model
    inits[[mean.model$fixed$name]] <- rep(0,ncol(jags.data$mean.fixed))
    if(!is.null(mean.model$random))
        inits[[mean.model$random$name]] <- rep(0,ncol(jags.data$mean.random))
    
    ## Generate initial values for parameters of variance model
    inits[[variance.model$fixed$name]] <- rep(0,ncol(jags.data$variance.fixed))
    if(!is.null(variance.model$random))
        inits[[variance.model$random$name]] <- rep(0,ncol(jags.data$variance.random))

    ## Initial response when rounding
    if(is.null(jags.data$response))
        inits$y <- (jags.data$lower + jags.data$upper)/2

    ## Return initial values list
    inits
}

    
