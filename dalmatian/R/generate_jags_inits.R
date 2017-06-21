generateJAGSinits <- function(mean.model,variance.model,jags.data,n.chains=1){

    ## Initialize list
    inits <- vector(mode="list",length=n.chains)

    for(i in 1:n.chains){
        ## Generate initial values for parameters of mean model
        inits[[i]][[mean.model$fixed$name]] <- rep(0,ncol(jags.data$mean.fixed))
        if(!is.null(mean.model$random))
            inits[[i]][[paste0(mean.model$random$name,".tmp")]] <- rep(0,ncol(jags.data$mean.random))

        ## Generate initial values for parameters of variance model
        inits[[i]][[variance.model$fixed$name]] <- rep(0,ncol(jags.data$variance.fixed))
        if(!is.null(variance.model$random))
            inits[[i]][[paste0(variance.model$random$name,".tmp")]] <- rep(0,ncol(jags.data$variance.random))

        ## Initial response when rounding
        if(is.null(jags.data$response))
            inits[[i]]$y <- (jags.data$lower + jags.data$upper)/2
    }

    ## Return initial values list
    inits
}

##' Set initial values for \code{dalmation}
##'
##' Allows the user to set initial values for \code{dalmation}. Any values
##' not specified will by initialized by \code{JAGS}.
##' @title Set initial values for \code{dalmation}
##' @param mean.model Model list specifying the structure of the mean. (list)
##' @param variance.model Model list specifyint the structure of the variance. (list)
##' @param fixed.mean Initial values for the fixed effects of the mean. (numeric)
##' @param fixed.variance Initial values for the fixed effects of the variance. (numeric)
##' @param y Inital values for the true response. This should only be specified if the \code{rounding = TRUE} in the main call to dalmation.
##' @param random.mean Initial values for the random effects of the mean. (numeric)
##' @param sd.mean Initial values for the standard deviation of the random effects of the mean. (numeric)
##' @param random.variance Initial values for the random effects of the variance. (numeric
##' @param sd.variance Initial values for the standard deviation of the random effects of the variance. (numeric)
##' @return inits (list)
##' @author Simon Bonner
##' @export
setJAGSInits <- function(mean.model,
                         variance.model,
                         fixed.mean=NULL,
                         fixed.variance=NULL,
                         y=NULL,
                         random.mean=NULL,
                         sd.mean=NULL,
                         random.variance=NULL,
                         sd.variance=NULL){

    ## Initialize list of initial values
    inits <- list()

    ## Set initial data values
    if(!is.null(y)){
      inits$y <- y
    }
    ## Set initial values for mean component of model
    ## 1) Fixed effects
    if(!is.null(fixed.mean)){
        inits[[mean.model$fixed$name]] <- fixed.mean
    }

    ## 2) Random effects
    if(!is.null(random.mean)){
        inits[[paste0(mean.model$random$name,".tmp")]] <- random.mean
    }

    ## 3) Random effects variances
    if(!is.null(sd.mean)){
        ## Generate redundant variables for Gelman's parametrization of half-t
        ncomp <- length(sd.mean)

        ## Compute variance parameter
        tau <- 1/sd.mean^2

        ## Set initial values
        inits[[paste0("redun.",mean.model$random$name)]] <- rep(1,ncomp)
        inits[[paste0("tau.",mean.model$random$name)]] <- tau
    }

    ## Set initial values for variance component of model
    ## 1) Fixed effects
    if(!is.null(fixed.variance)){
        inits[[variance.model$fixed$name]] <- fixed.variance
    }

    ## 2) Random effects
    if(!is.null(random.variance)){
        inits[[paste0(variance.model$random$name,".tmp")]] <- random.variance
    }

    ## 3) Random effects variance
    if(!is.null(sd.variance)){
        ## Generate redundant variables for Gelman's parametrization of half-t
        ncomp <- length(sd.variance)

        ## Compute variance parameter
        tau <- 1/sd.variance^2

        ## Set initial values
        inits[[paste0("redun.",variance.model$random$name)]] <- rep(1,ncomp)
        inits[[paste0("tau.",variance.model$random$name)]] <- tau
    }

    return(inits)
}


