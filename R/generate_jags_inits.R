##' @importFrom dglm dglm
##'
generateJAGSinits <- function(family,
                              mean.model,
                              dispersion.model,
                              jags.data){

  inits <- lapply(1:3,function(i){

    cat("     Initializing chain",i,"...\n")
    
    ## Initial response when rounding
    if(is.null(jags.data$y))
      y <- runif(jags.data$n,jags.data$lower,jags.data$upper)
    else
      y <- jags.data$y
    
    ## Remove zeros with negative binomial response and log link
    if(family == "nbinom" && mean.model$link == "log" && any(y == 0))
      y <- y + .1
    
    ## Mean formula
    if(is.null(mean.model$fixed) && is.null(mean.model$random)){
      stop("You have specified no fixed or random effects for the mean component of the model.\n\n")
    }
    else if(i==1 || is.null(mean.model$random)){ # Only fixed effects
      mean.formula <- formula("y ~ jags.data$mean.fixed - 1")
    }
    else if(i==2 || is.null(mean.model$fixed)){ # Only random effects
      mean.formula <- formula("y ~ jags.data$mean.random - 1")
    }
    else{ # Mixed effects
      mean.formula <- formula("y ~ jags.data$mean.fixed + jags.data$mean.random - 1")
    }
    
    # Dispersion formula
    if(is.null(dispersion.model$fixed) && is.null(dispersion.model$random)){
      stop("You have specified no fixed or random effects for the dispersion component of the model.\n\n")
    }
    else if(i==1 || is.null(dispersion.model$random)){ # Only fixed effects
      dispersion.formula <- formula("epsilonsq ~ jags.data$dispersion.fixed - 1")
    }
    else if(i==2 || is.null(dispersion.model$fixed)){ # Only random effects
      dispersion.formula <- formula("epsilonsq ~ jags.data$dispersion.random - 1")
    }
    else{ # Mixed effects
      dispersion.formula <- formula("epsilonsq ~ jags.data$dispersion.fixed + jags.data$dispersion.random - 1")
    }
    
    ##### My simple implementation of a double glm fit non-iteratively #####
    # ## Fit linear regression model
    # if(!is.null(mean.model$link))
    #   meanlm <- glm(mean.formula,family = gaussian(link=mean.model$link))
    # else
    #   meanlm <- lm(mean.formula)
    # 
    # ## Extract coefficients for fixed and random components of model
    # mean.coeff <- coef(meanlm)
    # 
    # fixed.mean <- mean.coeff[grep("fixed",names(mean.coeff))]
    # random.mean <- mean.coeff[grep("random",names(mean.coeff))]
    # 
    # ## Dispersion model
    # 
    # # Extract squared residuals from mean model
    # epsilonsq <- residuals(meanlm)^2
    # 
    # # Fit gamma GLM to squared residuals
    # dispersionlm <- glm(dispersion.formula,family=Gamma(link=dispersion.model$link))
    # 
    # # Extract coefficients
    # dispersion.coeff <- coef(dispersionlm)
    # fixed.dispersion <- dispersion.coeff[grep("fixed",names(dispersion.coeff))]
    # random.dispersion <- dispersion.coeff[grep("random",names(dispersion.coeff))]
    # 
    
    # Set link functions to identity if not specified
    if(is.null(mean.model$link))
      mean.model$link <- "identity"
    
    if(is.null(dispersion.model$link))
      dispersion.model$link <- "identity"
      
    # Fit double GLM (without random effects)
    dlink <- dispersion.model$link # I don't understand, but this is necessary.
    
    
    dglmfit <- dglm::dglm(formula=mean.formula,
                 dformula=dispersion.formula,
                 family=gaussian(link=mean.model$link),
                 dlink=dlink)
    
    # Extract coefficients of mean model
    mean.coeff <- coef(dglmfit)
 
    tmp <- grep("fixed",names(mean.coeff))
    if(length(tmp)>0)
      fixed.mean <- mean.coeff[tmp]
    else
      fixed.mean <- NULL
    
    tmp <- grep("random",names(mean.coeff))
    if(length(tmp)>0)
      random.mean <- mean.coeff[tmp]
    else
      random.mean <- NULL
    
    ## Compute random effects sd for mean
    if(!is.null(mean.model$random)){
      if(i %in% c(2,3)){
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
      else{
        ## Set random effects standard deviation to be very small
        ncomp <- jags.data[[paste0(mean.model$random$name,".ncomponents")]]
        
        sd.mean <- rep(.001,ncomp)
      }
    }
    else{
      sd.mean <- NULL
    }
    
    # Extract coefficients of dispersion model
    dispersion.coeff <- coef(dglmfit$dispersion)
    
    tmp <- grep("fixed",names(dispersion.coeff))
    if(length(tmp)>0)
      fixed.dispersion <- dispersion.coeff[tmp]
    else
      fixed.dispersion <- NULL
    
    tmp <- grep("random",names(dispersion.coeff))
    if(length(tmp)>0)
      random.dispersion <- dispersion.coeff[tmp]
    else
      random.dispersion <- NULL
    
    ## Compute random effects sd for dispersion
    if(!is.null(dispersion.model$random)){
      if(i %in% c(2,3)){
        ## Compute random effects standard deviations
        ncomp <- jags.data[[paste0(dispersion.model$random$name,".ncomponents")]]
        levels <- jags.data[[paste0(dispersion.model$random$name,".levels")]]
        
        sd.dispersion <- sapply(1:ncomp, function(j) sd(dispersion.coeff[which(levels==j)],na.rm=TRUE))
        
        ## Randomly fill in any missing random effects
        miss <- which(is.na(random.dispersion))
        
        if(length(miss) > 0){
          random.dispersion[miss] <- rnorm(length(miss),0,sd.dispersion[levels[miss]])
        }
      }
      else{
        ## Set random effects standard deviation to be very small
        ncomp <- jags.data[[paste0(dispersion.model$random$name,".ncomponents")]]
        
        sd.dispersion <- rep(.001,ncomp)
      }
    }
    else{
      sd.dispersion <- NULL
    }
    
    ## Construct initial values list
    if(is.null(jags.data$y))
      setJAGSInits(mean.model,
                   dispersion.model,
                   y=y,
                   fixed.mean=fixed.mean,
                   fixed.dispersion = fixed.dispersion,
                   random.mean=random.mean,
                   sd.mean=sd.mean,
                   random.dispersion=random.dispersion,
                   sd.dispersion = sd.dispersion)
    else
      setJAGSInits(mean.model,
                   dispersion.model,
                   fixed.mean=fixed.mean,
                   fixed.dispersion = fixed.dispersion,
                   random.mean=random.mean,
                   sd.mean=sd.mean,
                   random.dispersion=random.dispersion,
                   sd.dispersion = sd.dispersion)
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
##' @param dispersion.model Model list specifyint the structure of the dispersion. (list)
##' @param fixed.mean Initial values for the fixed effects of the mean. (numeric)
##' @param fixed.dispersion Initial values for the fixed effects of the dispersion. (numeric)
##' @param y Inital values for the true response. This should only be specified if the \code{rounding = TRUE} in the main call to dalmatian.
##' @param random.mean Initial values for the random effects of the mean. (numeric)
##' @param sd.mean Initial values for the standard deviation of the random effects of the mean. (numeric)
##' @param random.dispersion Initial values for the random effects of the dispersion. (numeric
##' @param sd.dispersion Initial values for the standard deviation of the random effects of the dispersion. (numeric)
##' @return inits (list)
##' @author Simon Bonner
##' @export
##' 
setJAGSInits <- function(mean.model,
                         dispersion.model,
                         fixed.mean=NULL,
                         fixed.dispersion=NULL,
                         y=NULL,
                         random.mean=NULL,
                         sd.mean=NULL,
                         random.dispersion=NULL,
                         sd.dispersion=NULL){

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

    ## 3) Random effects dispersions
    if(!is.null(sd.mean)){
        ## Generate redundant variables for Gelman's parametrization of half-t
        ncomp <- length(sd.mean)

        ## Compute dispersion parameter
        tau <- 1/sd.mean^2

        ## Set initial values
        inits[[paste0("redun.",mean.model$random$name)]] <- rep(1,ncomp)
        inits[[paste0("tau.",mean.model$random$name)]] <- tau
    }

    ## Set initial values for dispersion component of model
    ## 1) Fixed effects
    if(!is.null(fixed.dispersion)){
        inits[[dispersion.model$fixed$name]] <- fixed.dispersion
    }

    ## 2) Random effects
    if(!is.null(random.dispersion)){
        inits[[paste0(dispersion.model$random$name,".tmp")]] <- random.dispersion
    }

    ## 3) Random effects dispersion
    if(!is.null(sd.dispersion)){
        ## Generate redundant variables for Gelman's parametrization of half-t
        ncomp <- length(sd.dispersion)

        ## Compute dispersion parameter
        tau <- 1/sd.dispersion^2

        ## Set initial values
        inits[[paste0("redun.",dispersion.model$random$name)]] <- rep(1,ncomp)
        inits[[paste0("tau.",dispersion.model$random$name)]] <- tau
    }

    return(inits)
}


