generateJAGSlhd.gaussian <- function(jags.model.args,
                                     mean.model,
                                     dispersion.model,
                                     rounding,
                                     residuals,
                                     include.checks){
  
  ## Distribution of response
  string <- c("\t\t ## Data distribution \n",
              "\t\t y[i] ~ dnorm(muy[i],tauy[i])\n\n")

  ## Rounding 
  if(rounding){
    string <- c(string,
                generateJAGSrounding(jags.model.args))
  }

  ## Transform parameters
  string <- c(string,
              "\t\t ## Parameter transformations\n")
  
    if(is.null(dispersion.model$weights))
    string <- c(string,
                "\t\t tauy[i] <- 1/phi[i]\n")
  else
    string <- c(string,
                "\t\t tauy[i] <- weights[i]/phi[i]\n")

  ## Define variance and standard deviation as a function of mean and dispersion
  string <- c(string,
              "\t\t vary[i] <- phi[i]\n",
              "\t\t sdy[i] <- sqrt(vary[i])\n\n")

  ## Incorporate checks on the range of the mean and variance
  if(include.checks){
    string <- c(string,
                generateJAGSchecks(jags.model.args,
                                   "1", # Mean is any real value
                                   "(phi[i] > 0)")) # Variance must be positive
  }

  ## Pearson residuals
  if(residuals){
    string <- c(string,
                "\t\t ## Pearson residuals\n",
                "\t\t resid[i] <- (y[i] - muy[i])/sdy[i]\n\n")
  }

  string
}
  

