generateJAGSinits <- function(mean.model,variance.model,jags.data,n.chains=1){

    ## Initialize list
    inits <- vector(mode="list",length=n.chains)

    for(i in 1:n.chains){
        ## Generate initial values for parameters of mean model
        inits[[i]][[mean.model$fixed$name]] <- rep(0,ncol(jags.data$mean.fixed))
        if(!is.null(mean.model$random))
            inits[[i]][[mean.model$random$name]] <- rep(0,ncol(jags.data$mean.random))
        
        ## Generate initial values for parameters of variance model
        inits[[i]][[variance.model$fixed$name]] <- rep(0,ncol(jags.data$variance.fixed))
        if(!is.null(variance.model$random))
            inits[[i]][[variance.model$random$name]] <- rep(0,ncol(jags.data$variance.random))
        
        ## Initial response when rounding
        if(is.null(jags.data$response))
            inits[[i]]$y <- (jags.data$lower + jags.data$upper)/2
    }

    ## Return initial values list
    inits
}

    
