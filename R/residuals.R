##' Residuals method for dalmatian fitted objects
##'
##' Computes posterior summaries of the residuals for each observation.
##' Summary statistics include the posterior mean and the upper and
##' lower bounds of the 95% credible interval.
##' @title 
##' @param object Object of class \code{dalmatian} created by \code{dalmatian()}.
##' @param ... Ignored
##' @return Data frame containing original data augmented with posterior
##' mean and lower and upper bounds of the 95% credible interval of the
##' residual for each observation.
##' @export
##' @examples 
##' 
##' ## Load pied flycatcher data
##' data(pied_flycatchers_1)
##' 
##' ## Create variables bounding the true load
##' pfdata$lower=ifelse(pfdata$load==0,log(.001),log(pfdata$load-.049))
##' pfdata$upper=log(pfdata$load+.05)
##' 
##' ## Add 'log(IVI)' variable in pfdata
##' pfdata$'log(IVI)' <- log(pfdata$IVI)
##' 
##' ## Load output from previously run model
##' load(system.file("Pied_Flycatchers_1","pfresults.RData",package="dalmatian"))
##' 
##' ## Compute residuals
##' res.pfresults <- residual(object = pfresults)
##' 
##' @author Simon Bonner
residuals.dalmatian <- function(object,...){
  ## Extract residuals from MCMC sample (if computed)
  
  ## If residuals computed as part of MCMC simulation
  if(object$residuals){
    ## Extract sampled values of residuals
    cols <- grep("^resid\\[", colnames(object$coda[[1]]),value = TRUE)

    ## Compute summaries
    resid_summ <- summary(object$coda[,cols])

    ## Create output data frame
    output <- cbind(object$df,
                    Residual = resid_summ[[1]][,"Mean"],
                    Lower = resid_summ[[2]][,"2.5%"],
                    Upper = resid_summ[[2]][,"97.5%"])
  }
  ## Otherwise, if response is rounded
  else if(object$rounding){
    stop("If the response is rounded then the residuals must be computed ",
         "as part of the MCMC sampling. Please re-run the dalmatian() ",
         "command with residuals = TRUE.")
  }
  ## Otherwise, response is observed exactly
  else{
    ## Compute the fitted values for the mean
    output <- fitted(object)$mean

    ## Compute residuals
    y <- object$df[, object$response]

    output$Residual <- y - output$Fit
    output$Lower <- y - output$Lower
    otuput$Upper <- y - otuput$Upper
    output$Fit <- NULL
  }
  
  ## Return output
  output
}

  
