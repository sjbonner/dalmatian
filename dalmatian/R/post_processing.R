##' @importFrom stats start end window median sd
##' @importFrom coda thin
myCodaSummary <-
  function(coda,
           base,
           nstart = start(coda),
           nend = end(coda),
           nthin = coda::thin(coda)) {
    ## Generate summary for fixed effects components with names of the form base.xxx

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
          object$variance.model$fixed$name,
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

    if (!is.null(object$variance.model$random))
      output$varRandom = myCodaSummary(
        object$coda,
        paste0("sd\\.", object$variance.model$random$name),
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
  cat("Variance Model: Fixed Effects \n")
  print(round(x$varFixed, digits))

  if (!is.null(x$varRandom)) {
    cat("\n")
    cat("Variance Model: Random Effects \n")
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

    if (!is.null(object$variance.model$random))
      output$variance <-
        myCodaSummary(object$coda,
                      paste0("^", object$variance.model$random$name))

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
#' @param pars List of parameters to assess. If NULL (default) then diagnostics are computed for the fixed effects and random effects standard deviations in both the mean and variance models.
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
          paste0("^", object$variance.model$fixed$name, "\\."),
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

      if (!is.null(object$variance.model$random))
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
#' @param nstart Start point for computing summary statistics (relative to true start of chain).
#' @param nend End point for computing summary statistics (relative to true start of chain).
#' @param nthin Thinning factor for computing summary statsitics (relative to full chain and not previously thinned output).
#' @param family String defining selected family of variables (see help for \code{ggs()}).
#' @param plot If TRUE then generate plots. Otherwise, a list of \code{ggplot} objects will be returned.
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
           plot = TRUE,
           ...) {
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

      if (plot)
        print(output$meanFixed)

      ## Variance: fixed effects
      ggs2 <-
        ggmcmc::ggs(coda,
                    paste0("^", object$variance.model$fixed$name, "\\."))
      output$varianceFixed <- ggmcmc::ggs_traceplot(ggs2)

      if (plot)
        print(output$varianceFixed)

      ## Mean: random effects
      if (!is.null(object$mean.model$random)) {
        ggs3 <-
          ggmcmc::ggs(coda,
                      paste0("^sd\\.", object$mean.model$random$name))
        output$meanRandom <- ggmcmc::ggs_traceplot(ggs3)

        if (plot)
          print(output$meanRandom)
      }

      if (!is.null(object$variance.model$random)) {
        ggs4 <-
          ggmcmc::ggs(coda,
                      paste0("^sd\\.", object$variance.model$random$name))
        output$varianceRandom <- ggmcmc::ggs_traceplot(ggs4)

        if (plot)
          print(output$varianceRandom)
      }
    }
    else{
      ## Selected family
      ggs1 <- ggmcmc::ggs(object$coda, family)
      output <- ggmcmc::ggs_traceplot(ggs1)

      if (plot)
        print(output)
    }

    ## Return output
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
#' @param nstart Start point for computing summary statistics (relative to true start of chain).
#' @param nend End point for computing summary statistics (relative to true start of chain).
#' @param nthin Thinning factor for computing summary statsitics (relative to full chain and not previously thinned output).
#' @param family String defining selected family of variables (see help for \code{ggs()}).
#' @param plot If TRUE then generate plots. Otherwise, a list of \code{ggplot} objects will be returned.
#' @param ... Ignored
#'
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
           plot = TRUE,
           ...) {
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

      if (plot)
        print(output$meanFixed)

      ## Variance: fixed effects
      ggs2 <-
        ggmcmc::ggs(coda,
                    paste0("^", object$variance.model$fixed$name, "\\."))
      output$varianceFixed <- ggmcmc::ggs_caterpillar(ggs2)

      if (plot)
        print(output$varianceFixed)

      ## Mean: random effects
      if (!is.null(object$mean.model$random)) {
        ggs3 <-
          ggmcmc::ggs(coda,
                      paste0("^sd\\.", object$mean.model$random$name))
        output$meanRandom <- ggmcmc::ggs_caterpillar(ggs3)

        if (plot)
          print(output$meanRandom)
      }

      if (!is.null(object$variance.model$random)) {
        ggs4 <-
          ggmcmc::ggs(coda,
                      paste0("^sd\\.", object$variance.model$random$name))
        output$varianceRandom <- ggmcmc::ggs_caterpillar(ggs4)

        if (plot)
          print(output$varianceRandom)
      }
    }
    else{
      ## Selected family
      ggs1 <- ggmcmc::ggs(object$coda, family)
      output <- ggmcmc::ggs_caterpillar(ggs1)

      if (plot)
        print(output)
    }

    ## Return ouptut
    return(output)
  }
