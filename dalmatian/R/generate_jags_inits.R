generateJAGSinits <- function(mean.model,variance.model,jags.data,n.chains=1,spread=3){
  
  inits <- lapply(1:n.chains,function(i){
    ## Initial response when rounding
    if(is.null(jags.data$response))
      y <- runif(jags.data$n,jags.data$lower,jags.data$upper)
    else
      y <- jags.data$response
    
    ## Mean formula
    if(is.null(mean.model$fixed) && is.null(mean.model$random)){
      stop("You have specified no fixed or random effects for the mean component of the model.\n\n")
    }
    else if(is.null(mean.model$random)){ # Only fixed effects
      mean.formula <- formula("y ~ jags.data$mean.fixed - 1")
    }
    else if(is.null(mean.model$fixed)){ # Only random effects
      mean.formula <- formula("y ~ jags.data$mean.random - 1")
    }
    else{ # Mixed effects
      mean.formula <- formula("y ~ jags.data$mean.fixed + jags.data$mean.random - 1")
    }
    
    # Variance formula
    if(is.null(variance.model$fixed) && is.null(variance.model$random)){
      stop("You have specified no fixed or random effects for the random component of the model.\n\n")
    }
    else if(is.null(variance.model$random)){ # Only fixed effects
      variance.formula <- formula("epsilonsq ~ jags.data$variance.fixed - 1")
    }
    else if(is.null(variance.model$fixed)){ # Only random effects
      variance.formula <- formula("epsilonsq ~ jags.data$variance.random - 1")
    }
    else{ # Mixed effects
      variance.formula <- formula("epsilonsq ~ jags.data$variance.fixed + jags.data$variance.random - 1")
    }
    
    ##### My simple implementation of a double glm fit non-iteratively #####
    # ## Fit linear regression model
    # if(!is.null(mean.model$fixed$link))
    #   meanlm <- glm(mean.formula,family = gaussian(link=mean.model$fixed$link))
    # else
    #   meanlm <- lm(mean.formula)
    # 
    # ## Extract coefficients for fixed and random components of model
    # mean.coeff <- coef(meanlm)
    # 
    # fixed.mean <- mean.coeff[grep("fixed",names(mean.coeff))]
    # random.mean <- mean.coeff[grep("random",names(mean.coeff))]
    # 
    # ## Variance model
    # 
    # # Extract squared residuals from mean model
    # epsilonsq <- residuals(meanlm)^2
    # 
    # # Fit gamma GLM to squared residuals
    # variancelm <- glm(variance.formula,family=Gamma(link=variance.model$fixed$link))
    # 
    # # Extract coefficients
    # variance.coeff <- coef(variancelm)
    # fixed.variance <- variance.coeff[grep("fixed",names(variance.coeff))]
    # random.variance <- variance.coeff[grep("random",names(variance.coeff))]
    # 
    
    # Set link functions to identity if not specified
    if(is.null(mean.model$fixed$link))
      mean.model$fixed$link <- "identity"
    
    if(is.null(variance.model$fixed$link))
      variance.model$fixed$link <- "identity"
    
    # Fit double GLM (without random effects)
    dlink <- variance.model$fixed$link # I don't understand, but this is necessary.
    
    dglmfit <- dglm::dglm(formula=mean.formula,
                 dformula=variance.formula,
                 family=gaussian(link=mean.model$fixed$link),
                 dlink=dlink)
    
    # Extract coefficients of mean model
    mean.coeff <- coef(dglmfit)
    fixed.mean <- mean.coeff[grep("fixed",names(mean.coeff))]
    random.mean <- mean.coeff[grep("random",names(mean.coeff))]
    
    ## Compute random effects sd for mean
    if(!is.null(mean.model$random)){
      ## Compute random effects standard deviations
      ncomp <- jags.data[[paste0(mean.model$random$name,".ncomponents")]]
      levels <- jags.data[[paste0(mean.model$random$name,".levels")]]
      
      sd.mean <- sapply(1:ncomp, function(j) sd(mean.coeff[which(levels==j)],na.rm=TRUE))
      
      ## Randomly fill in any missing random effects
      miss <- which(is.na(random.mean))
      
      if(length(miss) > 0){
        random.mean[miss] <- rnorm(length(miss),0,sd.mean[levels[miss]])
      }
    }
    
    # Extract coefficients of variance model
    variance.coeff <- coef(dglmfit$dispersion)
    fixed.variance <- variance.coeff[grep("fixed",names(variance.coeff))]
    random.variance <- variance.coeff[grep("random",names(variance.coeff))]
    
    ## Compute random effects sd for variance
    if(!is.null(variance.model$random)){
      ## Compute random effects standard deviations
      ncomp <- jags.data[[paste0(variance.model$random$name,".ncomponents")]]
      levels <- jags.data[[paste0(variance.model$random$name,".levels")]]

      sd.variance <- sapply(1:ncomp, function(j) sd(variance.coeff[which(levels==j)],na.rm=TRUE))

      ## Randomly fill in any missing random effects
      miss <- which(is.na(random.variance))

      if(length(miss) > 0){
        random.variance[miss] <- rnorm(length(miss),0,sd.variance[levels[miss]])
      }
    }
    
    ## Construct initial values list
    if(is.null(jags.data$response))
      setJAGSInits(mean.model,
                   variance.model,
                   y=y,
                   fixed.mean=fixed.mean,
                   fixed.variance = fixed.variance,
                   random.mean=random.mean,
                   sd.mean=sd.mean,
                   random.variance=random.variance,
                   sd.variance = sd.variance)
    else
      setJAGSInits(mean.model,
                   variance.model,
                   fixed.mean=fixed.mean,
                   fixed.variance = fixed.variance,
                   random.mean=random.mean,
                   sd.mean=sd.mean,
                   random.variance=random.variance,
                   sd.variance = sd.variance)
  })
  
  ## Return initial values list
  inits
}

##' Set initial values for \code{dalmatian}
##'
##' Allows the user to set initial values for \code{dalmatian}. Any values
##' not specified will by initialized by \code{JAGS}.
##' @title Set initial values for \code{dalmatian}
##' @param mean.model Model list specifying the structure of the mean. (list)
##' @param variance.model Model list specifyint the structure of the variance. (list)
##' @param fixed.mean Initial values for the fixed effects of the mean. (numeric)
##' @param fixed.variance Initial values for the fixed effects of the variance. (numeric)
##' @param y Inital values for the true response. This should only be specified if the \code{rounding = TRUE} in the main call to dalmatian.
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


