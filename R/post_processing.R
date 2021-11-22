##' @importFrom coda thin
myCodaSummary <-
    function(coda,
             base,
             nstart = start(coda),
             nend = end(coda),
             nthin = coda::thin(coda)) {
        ## Generate summary for fixed effects components with names of the form base.
      
        ## Identify parameters matching the given form
        pars <-
            grep(paste0("^", base, "\\."), coda::varnames(coda), value = TRUE)

        ## Create new coda object with specified parameters
        ## Note: this is wasteful, but it is not possible to compute one set of HPD otherwise
        ## Note: It's also exactly what summary.mcmc.list does!
        coda1 <-
            coda::as.mcmc(do.call(rbind, window(
                                             coda[, pars, drop = FALSE],
                                             start = nstart,
                                             end = nend,
                                             thin = nthin
                                         )))

        ## Compute means and standard deviations
        basic <- cbind(apply(coda1, 2, mean),
                       apply(coda1, 2, median),
                       apply(coda1, 2, sd))

        ## Compute HPD intervals
        hpd1 <- coda::HPDinterval(coda1, prob = .95)
        hpd2 <- coda::HPDinterval(coda1, prob = .50)

        output <-
            cbind(basic, hpd1[, 1], hpd2[, 1], hpd2[, 2], hpd1[, 2])

        ## Add nice dimension names
        dimnames(output) <-
            list(
                sapply(pars, function(x)
                    strsplit(x, split = paste0("^", base, "."))[[1]][2]),
                c(
                    "Mean",
                    "Median",
                    "SD",
                    "Lower 95%",
                    "Lower 50%",
                    "Upper 50%",
                    "Upper 95%"
                )
            )

        ## Return object
        output
    }

#' Summary (dalmatian)
#'
#' @param object Object of class \code{dalmatian} created by \code{dalmatian()}.
#' @param nstart Start point for computing summary statistics (relative to true start of chain).
#' @param nend End point for computing summary statistics (relative to true start of chain).
#' @param nthin Thinning factor for computing summary statistics (relative to full chain and not previously thinned output).
#' @param ... Ignored
#'
#' @return output (list)
#' @references Bonner, S., Kim, H., Westneat, D., Mutzel, A.,
#' Wright, J., and Schofield, M.. (2021). \code{dalmatian}: A Package
#' for Fitting Double Hierarchical Linear Models in \code{R} via \code{JAGS} and
#' \code{nimble}. \emph{Journal of Statistical Software}, 100, 10, 1--25.
#' \doi{10.18637/jss.v100.i10}.
#' @export
#' 
summary.dalmatian <-
    function(object,
             nstart = start(object$coda),
             nend = end(object$coda),
             nthin = thin(object$coda),
             ...) {

      ## Initialize summary
      output <- list(start = nstart,
                     end = nend,
                     thin = nthin,
                     nchain = coda::nchain(object$coda))

      ## Compute summaries of fixed effects
      
      ## Mean
      if(!is.null(object$mean.model$fixed)){
        output$meanFixed <- myCodaSummary(
          object$coda,
          object$mean.model$fixed$name,
          nstart,
          nend,
          nthin
        )
      }

      ## Dispersion
      if(!is.null(object$dispersion.model$fixed)){
        output$dispFixed <- myCodaSummary(
          object$coda,
          object$dispersion.model$fixed$name,
          nstart,
          nend,
          nthin
        )
      }

      ## Joint
      if(!is.null(object$joint.model$fixed)){
        output$jointFixed <- myCodaSummary(
          object$coda,
          object$joint.model$fixed$name,
          nstart,
          nend,
          nthin
        )
      }
        
      
      ## Compute summaries of random effects

      ## Mean
        if (!is.null(object$mean.model$random))
            output$meanRandom = myCodaSummary(object$coda,
                                              paste0("sd\\.", object$mean.model$random$name),
                                              nstart,
                                              nend,
                                              nthin)

      ## Dispersion
        if (!is.null(object$dispersion.model$random))
            output$dispRandom = myCodaSummary(
                object$coda,
                paste0("sd\\.", object$dispersion.model$random$name),
                nstart,
                nend,
                nthin
            )

      ## Joint
      if (!is.null(object$joint.model$random))
            output$dispRandom = myCodaSummary(
                object$coda,
                paste0("sd\\.", object$joint.model$random$name),
                nstart,
                nend,
                nthin
            )

        class(output) <- "dalmatian.summary"

        return(output)
    }

#' Print Summary (dalmatian)
#'
#' @param x Object of class \code{dalamtion.summary} created by \code{summary.dalmatian()}.
#' @param digits Number of digits to display after decimal.
#' @param ... Ignored
#'
#' @return No return value. This function prints the summary of a dalmatian object in a nicely formatted manner. 
#' @references Bonner, S., Kim, H., Westneat, D., Mutzel, A.,
#' Wright, J., and Schofield, M.. (2021). \code{dalmatian}: A Package
#' for Fitting Double Hierarchical Linear Models in \code{R} via \code{JAGS} and
#' \code{nimble}. \emph{Journal of Statistical Software}, 100, 10, 1--25.
#' \doi{10.18637/jss.v100.i10}.
#' @export
#' 
print.dalmatian.summary <- function(x, digits = 2, ...) {
    ## Basic information about chains

    cat("\n", "Iterations = ", x$start, ":", x$end, "\n", sep = "")
    cat("Thinning interval =", x$thin, "\n")
    cat("Number of chains =", x$nchain, "\n")
    cat("Sample size per chain =",
    (x$end - x$start) / x$thin +
    1,
    "\n\n")

    cat("Posterior Summary Statistics for Each Model Component\n\n")

  if(!is.null(x$meanFixed)){
    cat("Mean Model: Fixed Effects \n")
    print(round(x$meanFixed, digits))
  }

  if (!is.null(x$meanRandom)) {
    cat("\n")
    cat("Mean Model: Random Effects \n")
    print(round(x$meanRandom, digits))
  }
  
  cat("\n")
  if (!is.null(x$dispFixed)){
    cat("Dispersion Model: Fixed Effects \n")
    print(round(x$dispFixed, digits))
  }

  if (!is.null(x$dispRandom)) {
    cat("\n")
    cat("Dispersion Model: Random Effects \n")
    print(round(x$dispRandom, digits))
  }
  
  cat("\n")
  if (!is.null(x$jointFixed)){
    cat("Joint Model: Fixed Effects \n")
    print(round(x$jointFixed, digits))
  }

  if (!is.null(x$jointRandom)) {
    cat("\n")
    cat("Joint Model: Random Effects \n")
    print(round(x$jointRandom, digits))
  }

  
}

#' Random Effects (S3 Generic)
#'
#' Generic function for exporting summaries of random effects.
#'
#' @param object Input object
#' @param ... Ignored
#'
#' @return List containing elements providing information on the predicted values of random effects as appropriate for the model. 
#' @references Bonner, S., Kim, H., Westneat, D., Mutzel, A.,
#' Wright, J., and Schofield, M.. (2021). \code{dalmatian}: A Package
#' for Fitting Double Hierarchical Linear Models in \code{R} via \code{JAGS} and
#' \code{nimble}. \emph{Journal of Statistical Software}, 100, 10, 1--25.
#' \doi{10.18637/jss.v100.i10}.
#' @export
#'
ranef <- function(object, ...) {
    UseMethod("ranef")
}

#' Random Effects (dalmatian)
#'
#' Compute posterior summary statistics for the individual random effects in each part of the model.
#'
#' @param object Object of class \code{dalmatian} created by \code{dalmatian()}.
#' @param nstart Start point for computing summary statistics (relative to true start of chain).
#' @param nend End point for computing summary statistics (relative to true start of chain).
#' @param nthin Thinning factor for computing summary statistics (relative to full chain and not previously thinned output).
#' @param ... Ignored
#'
#' @return List containing elements mean, dispersion, and/or joint as appropriate. Each element provides information on the predicted values of the random effects as appropriate for each component of the model. 
#' @references Bonner, S., Kim, H., Westneat, D., Mutzel, A.,
#' Wright, J., and Schofield, M.. (2021). \code{dalmatian}: A Package
#' for Fitting Double Hierarchical Linear Models in \code{R} via \code{JAGS} and
#' \code{nimble}. \emph{Journal of Statistical Software}, 100, 10, 1--25.
#' \doi{10.18637/jss.v100.i10}.
#' @export
#'
ranef.dalmatian <-
    function(object,
             nstart = start(object$coda),
             nend = end(object$coda),
             nthin = thin(object$coda),
             ...) {
        output <- list()

        if (!is.null(object$mean.model$random))
            output$mean <-
                myCodaSummary(object$coda, paste0("^", object$mean.model$random$name))

        if (!is.null(object$dispersion.model$random))
            output$dispersion <-
                myCodaSummary(object$coda,
                              paste0("^", object$dispersion.model$random$name))

        if (!is.null(object$joint.model$random))
          output$joint <-
            myCodaSummary(object$coda,
                          paste0("^", object$joint.model$random$name))

        return(output)
    }

#' Convergence Diagnostics (S3 Generic)
#'
#' Generic function for computing convergence diagnostics.
#'
#' @param object Object to asses.
#' @param ... Ignored
#'
#' @return List containing Gelman-Rubin and Raftery convergence diagnostics and effective sample sizes for the selected parameters. This information is used to diagnose convergence of the MCMC sampling algorithms.
#' @references Bonner, S., Kim, H., Westneat, D., Mutzel, A.,
#' Wright, J., and Schofield, M.. (2021). \code{dalmatian}: A Package
#' for Fitting Double Hierarchical Linear Models in \code{R} via \code{JAGS} and
#' \code{nimble}. \emph{Journal of Statistical Software}, 100, 10, 1--25.
#' \doi{10.18637/jss.v100.i10}.
#' @export
#'
convergence <- function(object, ...) {
    UseMethod("convergence")
}

#' Convergence
#'
#' Compute convergence diagnostics for a dalmatian object.
#'
#' @param object Object of class \code{dalmatian} created by \code{dalmatian()}.
#' @param pars List of parameters to assess. If NULL (default) then diagnostics are computed for the fixed effects and random effects standard deviations in the mean, dispersion, and joint components.
#' @param nstart Start point for computing summary statistics (relative to true start of chain).
#' @param nend End point for computing summary statistics (relative to true start of chain).
#' @param nthin Thinning factor for computing summary statistics (relative to full chain and not previously thinned output).
#' @param raftery List of arguments to be passed to \code{raftery.diag()}. Any values not provided will be set to their defaults (see \code{help(raftery.diag())} for details).
#' @param ... Ignored
#'
#' @return List containing Gelman-Rubin and Raftery convergence diagnostics and effective sample sizes for the selected parameters. This information is used to diagnose convergence of the MCMC sampling algorithms.
#' @references Bonner, S., Kim, H., Westneat, D., Mutzel, A.,
#' Wright, J., and Schofield, M.. (2021). \code{dalmatian}: A Package
#' for Fitting Double Hierarchical Linear Models in \code{R} via \code{JAGS} and
#' \code{nimble}. \emph{Journal of Statistical Software}, 100, 10, 1--25.
#' \doi{10.18637/jss.v100.i10}.
#' @export
#'
convergence.dalmatian <-
    function(object,
             pars = NULL,
             nstart = start(object$coda),
             nend = end(object$coda),
             nthin = coda::thin(object$coda),
             raftery=NULL,
             ...) {
        ## Select parameters to assess
        if (is.null(pars)) {
            pars <-
                c(grep(
                    paste0("^", object$mean.model$fixed$name, "\\."),
                    coda::varnames(object$coda),
                    value = TRUE
                ),
                grep(
                    paste0("^", object$dispersion.model$fixed$name, "\\."),
                    coda::varnames(object$coda),
                    value = TRUE
                ),
                grep(
                    paste0("^", object$joint.model$fixed$name, "\\."),
                    coda::varnames(object$coda),
                    value = TRUE
                ))

            if (!is.null(object$mean.model$random))
                pars <-
                    c(pars, grep(
                                paste0("sd\\.", object$mean.model$random$name),
                                coda::varnames(object$coda),
                                value = TRUE
                            ))

            if (!is.null(object$dispersion.model$random))
                pars <-
                    c(pars, grep(
                                paste0("sd\\.", object$dispersion.model$random$name),
                                coda::varnames(object$coda),
                                value = TRUE
                            ))

            if (!is.null(object$joint.model$random))
                pars <-
                    c(pars, grep(
                                paste0("sd\\.", object$joint.model$random$name),
                                coda::varnames(object$coda),
                                value = TRUE
                            ))
        }

        ## Set arguments for Raftery diagnostics
        if(is.null(raftery$q))
            raftery$q = .025
        if(is.null(raftery$r))
            raftery$r = .005
        if(is.null(raftery$s))
            raftery$s = .95
        if(is.null(raftery$converge.eps))
            raftery$converge.eps = .001

        ## Compute convergence diagnostics
        output <-
            list(gelman = coda::gelman.diag(window(
                                    object$coda[, pars],
                                    start = nstart,
                                    end = nend,
                                    thin = nthin,
                                    ),
                                    autoburnin = FALSE,
                                    multivariate = FALSE),
                 raftery = coda::raftery.diag(
                                     coda::as.mcmc(do.call(rbind, window(
                                                                      object$coda[, pars, drop = FALSE],
                                                                      start = nstart,
                                                                      end = nend,
                                                                      thin = nthin
                                                                  ))),
                                     q = raftery$q,
                                     r = raftery$r,
                                     s = raftery$s,
                                     converge.eps = raftery$converge.eps
                                 ),
                 effectiveSize=coda::effectiveSize(window(object$coda[,pars],
                                                          start=nstart,
                                                          end=nend,
                                                          thin=nthin)))

        return(output)
    }

#' Traceplots (Generic)
#'
#' @param object Object to assess.
#' @param ... Ignored
#'
#' @return A list of \code{ggplot} objects that can be used to later reproduce the plots via \code{print}.
#' @export
#' 
traceplots <- function(object, ...) {
    UseMethod("traceplots")
}

#' Traceplots (dalmatian)
#'
#' Construct traceplots for key (or selected) parameters in a dalmatian object.
#'
#' @param object Object of class \code{dalmatian} created by \code{dalmatian()}.
#' @param family String defining selected family of variables (see help for \code{ggs()}).
#' @param nstart Start point for computing summary statistics (relative to true start of chain).
#' @param nend End point for computing summary statistics (relative to true start of chain).
#' @param nthin Thinning factor for computing summary statistics (relative to full chain and not previously thinned output).
#' @param show If TRUE then plots are displayed on the computer screen and the session is paused between each plot. 
#' @param return_plots If TRUE then return list of \code{ggplot} objects.
#' @param ... Ignored
#'
#' @return A list of \code{ggplot} objects that can be used to later reproduce the plots via \code{print}.
#' @references Bonner, S., Kim, H., Westneat, D., Mutzel, A.,
#' Wright, J., and Schofield, M.. (2021). \code{dalmatian}: A Package
#' for Fitting Double Hierarchical Linear Models in \code{R} via \code{JAGS} and
#' \code{nimble}. \emph{Journal of Statistical Software}, 100, 10, 1--25.
#' \doi{10.18637/jss.v100.i10}.
#' @export
#'  
traceplots.dalmatian <-
    function(object,
             family = NULL,
             nstart = start(object$coda),
             nend = end(object$coda),
             nthin = thin(object$coda),
             show = TRUE,
             return_plots = TRUE,
             ...) {

        ## Identify whether session is interactive
        is_interactive <- interactive()

        ## If family is not specified then create all traceplots
        if (is.null(family)) {

            if(nstart != start(object$coda) || nend != end(object$coda) || nthin != thin(object$coda))
                coda <- window(object$coda,start=nstart,end=nend,thin=nthin)
            else
                coda <- object$coda
            
          ## Mean: fixed effects
          if(!is.null(object$mean.model$fixed)){
            ggs1 <-
              ggmcmc::ggs(coda,
                          paste0("^", object$mean.model$fixed$name, "\\."))
            output <- list(meanFixed = ggmcmc::ggs_traceplot(ggs1))

            if (show){
              if(is_interactive){
                readline(prompt="Press any key for the next plot:")
              }
              
              print(output$meanFixed)
            }
          }

          ## Dispersion: fixed effects
          if(!is.null(object$dispersion.model$fixed)){
            ggs2 <-
              ggmcmc::ggs(coda,
                          paste0("^", object$dispersion.model$fixed$name, "\\."))
            output$dispersionFixed <- ggmcmc::ggs_traceplot(ggs2)
            
            if (show){
              if(is_interactive){
                readline(prompt="Press any key for the next plot:")
              }
              
              print(output$dispersionFixed)
            }
          }
          
          ## Joint: fixed effects
          if(!is.null(object$joint.model$fixed)){
            ggs2 <-
              ggmcmc::ggs(coda,
                          paste0("^", object$joint.model$fixed$name, "\\."))
            output$jointFixed <- ggmcmc::ggs_traceplot(ggs2)
            
            if (show){
              if(is_interactive){
                readline(prompt="Press any key for the next plot:")
              }
              
              print(output$jointFixed)
            }
          }
          
          ## Mean: random effects
            if (!is.null(object$mean.model$random)) {
                ggs3 <-
                    ggmcmc::ggs(coda,
                                paste0("^sd\\.", object$mean.model$random$name))
                output$meanRandom <- ggmcmc::ggs_traceplot(ggs3)

                if (show){
                    if(is_interactive){
                        readline(prompt="Press any key for the next plot:")
                    }
                    
                    print(output$meanRandom)
                }
            }

          ## Dispersion: random effects
            if (!is.null(object$dispersion.model$random)) {
                ggs4 <-
                    ggmcmc::ggs(coda,
                                paste0("^sd\\.", object$dispersion.model$random$name))
                output$dispersionRandom <- ggmcmc::ggs_traceplot(ggs4)

                if (show){
                    if(is_interactive){
                        readline(prompt="Press any key for the next plot:")
                    }
                    
                    print(output$dispersionRandom)
                }
            }

          ## Joint: random effects
            if (!is.null(object$joint.model$random)) {
                ggs4 <-
                    ggmcmc::ggs(coda,
                                paste0("^sd\\.", object$joint.model$random$name))
                output$jointRandom <- ggmcmc::ggs_traceplot(ggs4)

                if (show){
                    if(is_interactive){
                        readline(prompt="Press any key for the next plot:")
                    }
                    
                    print(output$jointRandom)
                }
            }
        }
        else{
            ## Selected family
            ggs1 <- ggmcmc::ggs(object$coda, family)
            output <- ggmcmc::ggs_traceplot(ggs1)

            if (show){
                if(interactive()){
                    readline(prompt="Press any key for the next plot:")
                }
                
                print(output)
            }
        }

        ## Return output
        if(return_plots)
            output
    }

#' Caterpillar (Generic)
#'
#' @param object Object to assess.
#' @param ... Ignored
#'
#' @return A list of \code{ggplot} objects that can be used to later reproduce the plots via \code{print}.
#' @export
#'
caterpillar <- function(object, ...) {
    UseMethod("caterpillar")
}

#' Caterpillar (dalmatian)
#'
#' Construct caterpillar plots for key (or selected) parameters in a dalmatian object.
#'
#' @param object Object of class \code{dalmatian} created by \code{dalmatian()}.
#' @param family String defining selected family of variables (see help for \code{ggs()}).
#' @param nstart Start point for computing summary statistics (relative to true start of chain).
#' @param nend End point for computing summary statistics (relative to true start of chain).
#' @param nthin Thinning factor for computing summary statistics (relative to full chain and not previously thinned output).
#' @param show If TRUE then plots are displayed on the computer screen and the session is paused between each plot. 
#' @param return_plots If TRUE then a named list of ggplot objects containing the plots will be returned as output. 
#' @param ... Ignored
#'
#' @return A list of \code{ggplot} objects that can be used to later reproduce the plots via \code{print}.
#' @references Bonner, S., Kim, H., Westneat, D., Mutzel, A.,
#' Wright, J., and Schofield, M.. (2021). \code{dalmatian}: A Package
#' for Fitting Double Hierarchical Linear Models in \code{R} via \code{JAGS} and
#' \code{nimble}. \emph{Journal of Statistical Software}, 100, 10, 1--25.
#' \doi{10.18637/jss.v100.i10}.
#' @export
caterpillar.dalmatian <-
    function(object,
             family = NULL,
             nstart = start(object$coda),
             nend = end(object$coda),
             nthin = thin(object$coda),
             show = TRUE,
             return_plots = TRUE,
             ...) {

        ## Identify whether session is interactive
        is_interactive <- interactive()


        ## If family is not specified then create all caterpillar plots
        if (is.null(family)) {

            if(nstart != start(object$coda) || nend != end(object$coda) || nthin != thin(object$coda))
                coda <- window(object$coda,start=nstart,end=nend,thin=nthin)
            else
                coda <- object$coda

          ## Mean: fixed effects
          if(!is.null(object$mean.model$fixed)){
            ggs1 <-
              ggmcmc::ggs(coda,
                          paste0("^", object$mean.model$fixed$name, "\\."))
            output <- list(meanFixed = ggmcmc::ggs_caterpillar(ggs1))
            
            if (show){
              if(is_interactive){
                readline(prompt="Press any key for the next plot:")
              }
              
              print(output$meanFixed)
            }
          }

          ## Dispersion: fixed effects
          if(!is.null(object$dispersion.model$fixed)){
            ggs2 <-
              ggmcmc::ggs(coda,
                          paste0("^", object$dispersion.model$fixed$name, "\\."))
            output$dispersionFixed <- ggmcmc::ggs_caterpillar(ggs2)
            
            if (show){
              if(is_interactive){
                readline(prompt="Press any key for the next plot:")
              }
              
              print(output$dispersionFixed)
            }
          }
          
          ## Joint: fixed effects
          if(!is.null(object$joint.model$fixed)){
            ggs2 <-
              ggmcmc::ggs(coda,
                          paste0("^", object$joint.model$fixed$name, "\\."))
            output$jointFixed <- ggmcmc::ggs_caterpillar(ggs2)
            
            if (show){
              if(is_interactive){
                readline(prompt="Press any key for the next plot:")
              }
              
              print(output$jointFixed)
            }
          }
          
            ## Mean: random effects
          if (!is.null(object$mean.model$random)) {
                ggs3 <-
                    ggmcmc::ggs(coda,
                                paste0("^sd\\.", object$mean.model$random$name))
                output$meanRandom <- ggmcmc::ggs_caterpillar(ggs3)

                if (show){
                    if(is_interactive){
                        readline(prompt="Press any key for the next plot:")
                    }
                    
                    print(output$meanRandom)
                }
            }

          ## Dispersion: random effects
            if (!is.null(object$dispersion.model$random)) {
                ggs4 <-
                    ggmcmc::ggs(coda,
                                paste0("^sd\\.", object$dispersion.model$random$name))
                output$dispersionRandom <- ggmcmc::ggs_caterpillar(ggs4)

                if (show){
                    if(is_interactive){
                        readline(prompt="Press any key for the next plot:")
                    }
                    
                    print(output$dispersionRandom)
                }
            }

          ## Joint: random effects
            if (!is.null(object$joint.model$random)) {
                ggs4 <-
                    ggmcmc::ggs(coda,
                                paste0("^sd\\.", object$joint.model$random$name))
                output$jointRandom <- ggmcmc::ggs_caterpillar(ggs4)

                if (show){
                    if(is_interactive){
                        readline(prompt="Press any key for the next plot:")
                    }
                    
                    print(output$jointRandom)
                }
            }
        }
        else{
            ## Selected family
            ggs1 <- ggmcmc::ggs(object$coda, family)
            output <- ggmcmc::ggs_caterpillar(ggs1)

            if (show){
                if(is_interactive){
                    readline(prompt="Press any key for the next plot:")
                }
                
                print(output)
            }
        }

        ## Return output
        if(return_plots)
            output

    }
##' terms (dalmatian)
##'
##' Constructs a list of terms objects for each component of the model specified in the
##' input object.
##' 
##' @title Terms function for \code{dalmatian} objects
##' @param x Object of class \code{dalmatian} created by \code{dalmatian()}.
##' @param ... Further object passed directly to \code{terms}. Recycled for each model component.
##' @return List of with two lists named mean and dispersion each containing \code{terms} objects
##' corresponding to the fixed and random components of that model component (if present).
##' @export
##' @author Simon Bonner
terms.dalmatian <- function(x,...){
  ## Extract terms objects for each component of the fitted model

  terms.dalmatian.component <- function(model,...){
    ## Local function to extract terms form mean or dispersion component individually

    ## a) Fixed effects
    if (!is.null(model$fixed))
      terms_fixed <- terms(model$fixed$formula, ...)
    else
      terms_fixed <- NULL

    ## b) Random effects
    if (!is.null(model$random))
      terms_random <- terms(model$random$formula, ...)
    else
      terms_random <- NULL

    ## Return
    list(fixed = terms_fixed,
         random = terms_random)
  }
  
  ## Return output
  list(mean = terms.dalmatian.component(x$mean.model,...),
       dispersion = terms.dalmatian.component(x$dispersion.model,...),
       joint = terms.dalmatian.component(x$joint.model,...))
}

##' coef (dalmatian)
##'
##' Extracts coefficients for the mean and dispersion components of a
##' dalmatian model.
##' 
##' @title Coefficients function for \code{dalmatian} objects
##' @param object Object of class \code{dalmatian} created by \code{dalmatian()}.
##' @param summary Posterior summaries computed from the supplied \code{dalmatian} object (optional).
##' @param ranef Random effects summary computed from the supplied \code{dalmatian} object (optional).
##' @param ... Ignored
##' @return List of three lists named mean, dispersion, and joint each containing the posterior means of the coefficients
##' corresponding to the fixed and random terms of that model component (if present).
##' @references Bonner, S., Kim, H., Westneat, D., Mutzel, A.,
##' Wright, J., and Schofield, M.. (2021). \code{dalmatian}: A Package
##' for Fitting Double Hierarchical Linear Models in \code{R} via \code{JAGS} and
##' \code{nimble}. \emph{Journal of Statistical Software}, 100, 10, 1--25.
##' \doi{10.18637/jss.v100.i10}.
##' @export
##' @author Simon Bonner
coef.dalmatian <- function(object,summary = NULL, ranef = NULL, ...){
  
  ## Compute posterior summaries if not provided'
  if(is.null(summary))
    summary <- summary(object)

  ## Compute posterior summaries of random effects (if not provided and model contains random effects)
  if((!is.null(object$mean.model$random) | !is.null(object$dispersion.model$random) |
      !is.null(object$joint.model$random)) & is.null(ranef))
    ranef <- ranef(object)

  ## Mean model
  
  ## Extract posterior means for fixed effects
  mean_fixef <- summary$meanFixed[,"Mean"]

  ## If model only contains fixed effects
  if(is.null(object$mean.model$random))
    coef_mean <- mean_fixef
  
  ## Otherwise, combine fixed and random effects
  else{
    ## Extract and format posterior means for random effects
    mean_ranef <- dplyr::as_tibble(ranef$mean,rownames = "Effect") %>%
      dplyr::select(.data$Effect, .data$Mean) %>%
      tidyr::separate(.data$Effect, c("ID","Effect"),sep=":",fill="right") %>%
      tidyr::replace_na(list(Effect = "(Intercept)")) %>%
      tidyr::spread(key = .data$Effect, value = .data$Mean)
    
    ## Combine fixed and random effects
    allef <- unique(c(names(mean_fixef),names(mean_ranef)[-1]))
    
    tmp <- lapply(allef,function(ef){
      if(!ef %in% names(mean_fixef))
        dplyr::pull(mean_ranef,ef)
      else if(!ef %in% names(mean_ranef))
        rep(mean_fixef[ef],nrow(mean_ranef))
      else
        mean_fixef[ef] + dplyr::pull(mean_ranef,ef)
    }) 

    coef_mean <- do.call("cbind",tmp)

    ## Add appropriate dimension names
    dimnames(coef_mean) <- list(dplyr::pull(mean_ranef,"ID"),allef)

    ## Return output
    coef_mean
  }

  ## Dispersion model
  
  ## Extract posterior means for fixed effects
  disp_fixef <- summary$dispFixed[,"Mean"]

  ## If model only contains fixed effects
  if(is.null(object$dispersion.model$random))
    coef_disp <- disp_fixef
  
  ## Otherwise, combine fixed and random effects
  else{
    ## Extract and format posterior means for random effects
    disp_ranef <- dplyr::as_tibble(ranef$disp,rownames = "Effect") %>%
      dplyr::select(.data$Effect, .data$Mean) %>%
      tidyr::separate(.data$Effect, c("ID","Effect"),sep=":",fill="right") %>%
      tidyr::replace_na(list(Effect = "(Intercept)")) %>%
      tidyr::spread(key = .data$Effect, value = .data$Mean)
    
    ## Combine fixed and random effects
    allef <- unique(c(names(disp_fixef),names(disp_ranef)[-1]))
    
    tmp <- lapply(allef,function(ef){
      if(!ef %in% names(disp_fixef))
        dplyr::pull(disp_ranef,ef)
      else if(!ef %in% names(disp_ranef))
        rep(disp_fixef[ef],nrow(disp_ranef))
      else
        disp_fixef[ef] + dplyr::pull(disp_ranef,ef)
    }) 

    coef_disp <- do.call("cbind",tmp)

    ## Add appropriate dimension names
    dimnames(coef_disp) <- list(dplyr::pull(disp_ranef,"ID"),allef)

    ## Return output
    coef_disp
  }

  ## Joint model
  
  ## Extract posterior means for fixed effects
  joint_fixef <- summary$jointFixed[,"Mean"]

  ## If model only contains fixed effects
  if(is.null(object$joint.model$random))
    coef_joint <- joint_fixef
  
  ## Otherwise, combine fixed and random effects
  else{
    ## Extract and format posterior means for random effects
    joint_ranef <- dplyr::as_tibble(ranef$joint,rownames = "Effect") %>%
      dplyr::select(.data$Effect, .data$Mean) %>%
      tidyr::separate(.data$Effect, c("ID","Effect"),sep=":",fill="right") %>%
      tidyr::replace_na(list(Effect = "(Intercept)")) %>%
      tidyr::spread(key = .data$Effect, value = .data$Mean)
    
    ## Combine fixed and random effects
    allef <- unique(c(names(joint_fixef),names(joint_ranef)[-1]))
    
    tmp <- lapply(allef,function(ef){
      if(!ef %in% names(joint_fixef))
        dplyr::pull(joint_ranef,ef)
      else if(!ef %in% names(joint_ranef))
        rep(joint_fixef[ef],nrow(joint_ranef))
      else
        joint_fixef[ef] + dplyr::pull(joint_ranef,ef)
    }) 

    coef_joint <- do.call("cbind",tmp)

    ## Add appropriate dimension names
    dimnames(coef_joint) <- list(dplyr::pull(joint_ranef,"ID"),allef)

    ## Return output
    coef_joint
  }
  
  ## Return output
  list(mean = coef_mean,
       dispersion = coef_disp,
       joint = coef_joint)
}
