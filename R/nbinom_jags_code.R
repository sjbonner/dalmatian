generateJAGSlhd.nbinom <- function(jags.model.args,
                                   mean.model,
                                   dispersion.model,
                                   rounding,
                                   residuals,
                                   include.checks){
  
  ## Distribution of response
  string <- c("\t\t ## Data distribution \n",
              "\t\t y[i] ~ dnegbin(p[i],r[i])\n\n")

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
                "\t\t r[i] <- 1/phi[i]\n",
                "\t\t p[i] <- 1/(1 + phi[i] * muy[i])\n")
  else
    stop("Weights are currently not supported for the negative binomial model.\n")

  ## Define variance and standard deviation as a function of mean and dispersion
  string <- c(string,
              "\t\t vary[i] <- muy[i] / p[i]\n",
              "\t\t sdy[i] <- sqrt(vary[i])\n\n")

  ## Incorporate checks on the range of the mean and variance
  if(include.checks){
    string <- c(string,
                generateJAGSchecks(jags.model.args,
                                   "(muy[i] > 0)", # Mean must be positive
                                   "(phi[i] > 0)")) # Dispersion must be positive
  }

  ## Pearson residuals
  if(residuals){
    string <- c(string,
                "\t\t ## Pearson residuals\n",
                "\t\t resid[i] <- (y[i] - muy[i])/sdy[i]\n\n")
  }

  string
}
  
