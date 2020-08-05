##' @importFrom stats start end window median sd
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
#' @param nthin Thinning factor for computing summary statsitics (relative to full chain and not previously thinned output).
#' @param ... Ignored
#'
#' @return output (list)
#' @export
#' 
#' @examples 
#' 
#' ## Load output from previously run model
#' load(system.file("Pied_Flycatchers_1","pfresults.RData",package="dalmatian"))
#' 
#' ## Compute numerical summaries
#' summary(pfresults)
summary.dalmatian <-
    function(object,
             nstart = start(object$coda),
             nend = end(object$coda),
             nthin = thin(object$coda),
             ...) {
        ## Compute summaries of fixed effects
        output <-
            list(
                meanFixed = myCodaSummary(
                    object$coda,
                    object$mean.model$fixed$name,
                    nstart,
                    nend,
                    nthin
                ),
                varFixed = myCodaSummary(
                    object$coda,
                    object$dispersion.model$fixed$name,
                    nstart,
                    nend,
                    nthin
                ),
                start = nstart,
                end = nend,
                thin = nthin,
                nchain = coda::nchain(object$coda)
            )

        ## Compute summaries of random effects
        if (!is.null(object$mean.model$random))
            output$meanRandom = myCodaSummary(object$coda,
                                              paste0("sd\\.", object$mean.model$random$name),
                                              nstart,
                                              nend,
                                              nthin)

        if (!is.null(object$dispersion.model$random))
            output$varRandom = myCodaSummary(
                object$coda,
                paste0("sd\\.", object$dispersion.model$random$name),
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
#' @export
#' 
#' @examples 
#' 
#' ## Load output from previously run model
#' load(system.file("Pied_Flycatchers_1","pfresults.RData",package="dalmatian"))
#' 
#' ## Compute numerical summaries
#' print(summary(pfresults))
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

    cat("Mean Model: Fixed Effects \n")
    print(round(x$meanFixed, digits))

    if (!is.null(x$meanRandom)) {
        cat("\n")
        cat("Mean Model: Random Effects \n")
        print(round(x$meanRandom, digits))
    }

    cat("\n")
    cat("Dispersion Model: Fixed Effects \n")
    print(round(x$varFixed, digits))

    if (!is.null(x$varRandom)) {
        cat("\n")
        cat("Dispersion Model: Random Effects \n")
        print(round(x$varRandom, digits))
    }
}

#' Random Effects (S3 Generic)
#'
#' Generic function for exporting summaries of random effects.
#'
#' @param object Input object
#' @param ... Ignored
#'
#' @export
#'
#' @examples 
#' 
#' ## Load output from previously run model
#' load(system.file("Pied_Flycatchers_1","pfresults.RData",package="dalmatian"))
#' 
#' ## Compute numerical summaries
#' ranef(pfresults)
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
#' @param nthin Thinning factor for computing summary statsitics (relative to full chain and not previously thinned output).
#' @param ... Ignored
#'
#' @return output (list)
#' @export
#'
#' @examples 
#' 
#' ## Load output from previously run model
#' load(system.file("Pied_Flycatchers_1","pfresults.RData",package="dalmatian"))
#' 
#' ## Compute numerical summaries
#' ranef(pfresults)
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

        return(output)
    }

#' Convergence Diagnostics (S3 Generic)
#'
#' Generic function for computing convergence diagnostics.
#'
#' @param object Object to asses.
#' @param ... Ignored
#'
#' @export
#'
#' @examples 
#' 
#' ## Load output from previously run model
#' load(system.file("Pied_Flycatchers_1","pfresults.RData",package="dalmatian"))
#' 
#' ## Compute convergence diagnostics
#' pfconvergence <- convergence(pfresults)
#'
convergence <- function(object, ...) {
    UseMethod("convergence")
}

#' Convergence
#'
#' Compute convergence diagnostics for a dalmatian object.
#'
#' @param object Object of class \code{dalmatian} created by \code{dalmatian()}.
#' @param pars List of parameters to assess. If NULL (default) then diagnostics are computed for the fixed effects and random effects standard deviations in both the mean and dispersion models.
#' @param nstart Start point for computing summary statistics (relative to true start of chain).
#' @param nend End point for computing summary statistics (relative to true start of chain).
#' @param nthin Thinning factor for computing summary statsitics (relative to full chain and not previously thinned output).
#' @param raftery List of arguments to be passed to \code{raftery.diag()}. Any values not provided will be set to their defaults (see \code{help(raftery.diag())} for details).
#' @param ... Ignored
#'
#' @return List containing Gelman-Rubin and Raftery convergence diagnostics and effective sampel sizes for the selected parameters.
#' @export
#'
#' @examples 
#' 
#' ## Load output from previously run model
#' load(system.file("Pied_Flycatchers_1","pfresults.RData",package="dalmatian"))
#' 
#' ## Compute convergence diagnostics
#' pfconvergence <- convergence(pfresults)
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
                                paste0("sd\\.", object$variace.model$random$name),
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
#' @export
#' 
#' @examples 
#' 
#' ## Load output from previously run model
#' load(system.file("Pied_Flycatchers_1","pfresults.RData",package="dalmatian"))
#' 
#' ## Generate traceplots
#' pftraceplots <- traceplots(pfresults)
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
#' @param nthin Thinning factor for computing summary statsitics (relative to full chain and not previously thinned output).
#' @param show If TRUE then plots are displayed on the computer screen and the session is paused between each plot. 
#' @param return_plots If TRUE then return list of \code{ggplot} objects.
#' @param ... Ignored
#'
#' @return A list of \code{ggplot} objects that can be used to later reproduce the plots via \code{print}.
#' @export
#'  
#' @examples 
#' 
#' ## Load output from previously run model
#' load(system.file("Pied_Flycatchers_1","pfresults.RData",package="dalmatian"))
#' 
#' ## Generate traceplots
#' pftraceplots <- traceplots(pfresults)
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

            ## Dispersion: fixed effects
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
#' @export
#'
#' @examples 
#' 
#' ## Load output from previously run model
#' load(system.file("Pied_Flycatchers_1","pfresults.RData",package="dalmatian"))
#' 
#' ## Generate caterpillar
#' pfcaterpillar <- caterpillar(pfresults,plot = FALSE)
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
#' @param nthin Thinning factor for computing summary statsitics (relative to full chain and not previously thinned output).
#' @param show If TRUE then plots are displayed on the computer screen and the session is paused between each plot. 
#' @param return_plots If TRUE then a named list of ggplot objects containing the plots will be returned as output. 
#' @param ... Ignored
#'
#' @param plot 
#' @return A list of \code{ggplot} objects that can be used to later reproduce the plots via \code{print}.
#' @export
#' @examples 
#' 
#' ## Load output from previously run model
#' load(system.file("Pied_Flycatchers_1","pfresults.RData",package="dalmatian"))
#' 
#' ## Generate caterpillar
#' pfcaterpillar <- caterpillar(pfresults,plot = FALSE)
#'
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

            ## Dispersion: fixed effects
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
##' @examples
##' \dontrun{
##' ## Extract the terms objects corresponding to the pied-flycatcher model without random effects
##' terms(pfresults)
##'
##' ## Extract the terms objects corresponding to the pied-flycatcher model with random effects
##' terms(pfresults2)
##' }
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
       dispersion = terms.dalmatian.component(x$dispersion.model,...))
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
##' @return List of two lists named mean and dispersion each containing the posterior means of the coefficients
##' corresponding to the fixed and random terms of that model component (if present).
##' @author Simon Bonner
coef.dalmatian <- function(object,summary = NULL, ranef = NULL){
  ## Compute posterior summaries if not provided'
  if(is.null(summary))
    summary <- summary(object)

  ## Compute posterior summaries of random effects (if not provided and model contains random effects)
  if((!is.null(object$mean.model$random) | !is.null(object$dispersion.model$random)) & is.null(ranef))
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
  var_fixef <- summary$varFixed[,"Mean"]

  ## If model only contains fixed effects
  if(is.null(object$var.model$random))
    coef_var <- var_fixef
  
  ## Otherwise, combine fixed and random effects
  else{
    ## Extract and format posterior means for random effects
    var_ranef <- dplyr::as_tibble(ranef$var,rownames = "Effect") %>%
      dplyr::select(.data$Effect, .data$Mean) %>%
      tidyr::separate(.data$Effect, c("ID","Effect"),sep=":",fill="right") %>%
      tidyr::replace_na(list(Effect = "(Intercept)")) %>%
      tidyr::spread(key = .data$Effect, value = .data$Mean)
    
    ## Combine fixed and random effects
    allef <- unique(c(names(var_fixef),names(var_ranef)[-1]))
    
    tmp <- lapply(allef,function(ef){
      if(!ef %in% names(var_fixef))
        dplyr::pull(var_ranef,ef)
      else if(!ef %in% names(var_ranef))
        rep(var_fixef[ef],nrow(var_ranef))
      else
        var_fixef[ef] + dplyr::pull(var_ranef,ef)
    }) 

    coef_var <- do.call("cbind",tmp)

    ## Add appropriate dimension names
    dimnames(coef_var) <- list(dplyr::pull(var_ranef,"ID"),allef)

    ## Return output
    coef_var
  }

  ## Return output
  list(mean = coef_mean,
       dispersion = coef_var)
}
