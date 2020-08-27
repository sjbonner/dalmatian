generateJAGSlhd.betabin <- function(jags.model.args,
                                    mean.model,
                                    dispersion.model,
                                    rounding,
                                    residuals,
                                    include.checks){
  
  ## Distribution of response
  string <- c("\t\t ## Data distribution \n",
              "\t\t y[i] ~ dbetabin(alphay[i], betay[i], m[i])\n\n")

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
                "\t\t alphay[i] <- muy[i] * (1 - phi[i]) / phi[i]\n",
                "\t\t betay[i] <- (1 - muy[i]) * (1 - phi[i]) / phi[i]\n")
  else
    stop("Weights are currently not supported for the beta-binomial model.\n")

  ## Define variance and standard deviation as a function of mean and dispersion
  string <- c(string,
              "\t\t vary[i] <- m[i] * muy[i] * (1-muy[i]) * (1 + (m[i] - 1) * phi[i])\n",
              "\t\t sdy[i] <- sqrt(vary[i])\n\n")

  ## Incorporate checks on the range of the mean and variance
  if(include.checks){
    string <- c(string,
                generateJAGSchecks(jags.model.args,
                                   "(muy[i] > 0) * (muy[i] < 1)", # Mean must be between 0 and 1
                                   "(phi[i] > 0) * (phi[i] < 1)")) # Dispersion must be between 0 and 1
  }

  ## Pearson residuals
  if(residuals){
    string <- c(string,
                "\t\t ## Pearson residuals\n",
                "\t\t resid[i] <- (y[i] - m[i] * muy[i])/sdy[i]\n\n")
  }

  string
}
  
