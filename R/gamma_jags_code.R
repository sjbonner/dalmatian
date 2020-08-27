generateJAGSlhd.gamma <- function(jags.model.args,
                                  mean.model,
                                  dispersion.model,
                                  rounding,
                                  residuals,
                                  include.checks){
  
  ## Distribution of response
  string <- c("\t\t ## Data distribution \n",
              "\t\t y[i] ~ dgamma(ry[i], lambday[i])\n\n")

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
                "\t\t lambday[i] <- weights[i]/phi[i]\n",
                "\t\t ry[i] <- weights[i] * muy[i] * lambday[i]\n")
  else
    string <- c(string,
                "\t\t lambday[i] <- 1/phi[i]\n",
                "\t\t ry[i] <- muy[i] * lambday[i]\n")

  ## Define variance and standard deviation as a function of mean and dispersion
  string <- c(string,
              "\t\t vary[i] <- phi[i] * muy[i]\n",
              "\t\t sdy[i] <- sqrt(vary[i])\n\n")

  ## Incorporate checks on the range of the mean and variance
  if(include.checks){
    string <- c(string,
                generateJAGSchecks(jags.model.args,
                                   "(muy[i] > 0)", # Mean must be positive
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
  
