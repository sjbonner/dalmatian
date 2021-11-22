##' The primary function which automates the running of \code{JAGS} and \code{nimble}.
##'
##' The primary function in the package, dalmatian automates the generation of code, data, and initial values. These are then passed as arguments to function from the \code{rjags} package which automates the generation of samples from the posterior.
##' @title Run DGLM in \code{JAGS} via \code{rjags} or in \code{nimble}
##'
##' @param df Data frame containing the response and predictor values for each individual. (data.frame)
##' @param family Name of family of response distribution. Currently supported families include normal (\code{gaussian}) and negative binomial (\code{nbinom}). (character)
##' @param mean.model Model list specifying the structure of the mean. (list)
##' @param dispersion.model Model list specifying the structure of the dispersion. (list)
##' @param joint.model Model list specifying structure with parameter shared between linear predictors of the mean and variance. (list)
##' @param jags.model.args List containing named arguments of \code{jags.model}. (list)
##' @param coda.samples.args List containing named arguments of \code{coda.samples}. (list)
##' @param response Name of variable in the data frame representing the response. (character)
##' @param ntrials Name of variable in the data frame representing the number of independent trials for each observation of the beta binomial model.
##' @param rounding Specifies that response has been rounded if TRUE. (logical)
##' @param lower Name of variable in the data frame representing the lower bound on the response if rounded. (character)
##' @param upper Name of variable in the data frame representing the upper bound on the response if rounded. (character)
##' @param parameters Names of parameters to monitor. If NULL then default values are selected. (character)
##' @param svd Compute Singular Variable Decomposition of model matrices to improve convergence. (logical)
##' @param residuals If TRUE then compute residuals in output. (logical)
##' @param gencode If TRUE then generate code potentially overwriting existing model file. By default generate code if the file does not exist and prompt user if it does. (logical)
##' @param run.model If TRUE then run sampler. Otherwise, stop once code and data have been created. (logical)
##' @param engine Specifies the sampling software. Packages currently supported include JAGS (the default) and nimble. (character)
##' @param n.cores Number of cores to use. If equal to 1 then chains will not be run in parallel. If greater than 1 then chains will be run in parallel using the designated number of cores.
##' @param include.checks If TRUE (default) then include extra Bernoulli variables in the model to ensure that the mean and dispersion parameters remain within their support. (logical)
##' @param drop.levels If TRUE then drop unused levels from all factors in \code{df}. (logical)
##' @param drop.missing If TRUE then remove records with missing response variable. (logical)
##' @param overwrite If TRUE then overwrite existing JAGS files (non-interactive sessions only). (logical)
##' @param debug If TRUE then enter debug model. (logical)
##' @param saveJAGSinput Directory to which jags.model input is saved prior to calling \code{jags.model()}. This is useful for debugging. No files saved if NULL. (character)
##'
##' @return An object of class \code{dalmatian} containing copies of the original data frame, the mean model, the
##' dispersion model the arguments of \code{jags.model} and \code{coda.samples}. and the output of the MCMC sampler. 
##' @author Simon Bonner
##' @references Bonner, S., Kim, H., Westneat, D., Mutzel, A.,
##' Wright, J., and Schofield, M.. (2021). \code{dalmatian}: A Package
##' for Fitting Double Hierarchical Linear Models in \code{R} via \code{JAGS} and
##' \code{nimble}. \emph{Journal of Statistical Software}, 100, 10, 1--25.
##' \doi{10.18637/jss.v100.i10}.
##' @export
##' @examples
##' 
##' \dontrun{
##' ## Load pied flycatcher data
##' data(pied_flycatchers_1)
##' 
##' ## Create variables bounding the true load
##' pfdata$lower=ifelse(pfdata$load==0,log(.001),log(pfdata$load-.049))
##' pfdata$upper=log(pfdata$load+.05)

##' ## Mean model
##' mymean=list(fixed=list(name="alpha",
##'                       formula=~ log(IVI) + broodsize + sex,
##'                       priors=list(c("dnorm",0,.001))))
##'
##' ## Dispersion model
##' myvar=list(fixed=list(name="psi",
##'                       link="log",
##'                       formula=~broodsize + sex,
##'                       priors=list(c("dnorm",0,.001))))
##' 
##' ## Set working directory
##' ## By default uses a system temp directory. You probably want to change this.
##' workingDir <- tempdir()
##' 
##' ## Define list of arguments for jags.model()
##' jm.args <- list(file=file.path(workingDir,"pied_flycatcher_1_jags.R"),n.adapt=1000)
##' 
##' ## Define list of arguments for coda.samples()
##' cs.args <- list(n.iter=5000)
##' 
##' ## Run the model using dalmatian
##' pfresults <- dalmatian(df=pfdata,
##'                          mean.model=mymean,
##'                          dispersion.model=myvar,
##'                          jags.model.args=jm.args,
##'                          coda.samples.args=cs.args,
##'                          rounding=TRUE,
##'                          lower="lower",
##'                          upper="upper",
##'                          debug=FALSE)
##' }                          
dalmatian <- function(df,
                      family = "gaussian",
                      mean.model,
                      dispersion.model,
                      joint.model = NULL,
                      jags.model.args,
                      coda.samples.args,
                      response = NULL,
                      ntrials = NULL,
                      rounding = FALSE,
                      lower = NULL,
                      upper = NULL,
                      parameters = NULL,
                      svd = FALSE,
                      residuals = FALSE,
                      gencode = NULL,
                      run.model = TRUE,
                      engine = "JAGS",
                      n.cores = 1L,
                      drop.levels = TRUE,
                      drop.missing = TRUE,
                      include.checks = TRUE,
                      overwrite = FALSE,
                      debug = FALSE,
                      saveJAGSinput=NULL) {
  
  ## Enter debug state
  if (debug)
    browser()

  ## Check that input is consistent and sufficient
  if (! family %in% c("gaussian","nbinom","betabin","gamma"))
    stop("Currently supported families of distributions for the response include either Gaussian (family=\"gaussian\"), negative binomial (family=\"nbinom\"), beta-binomial (family = \"betabin\") of gamma (family = \"gamma\").\n\n")
  
  if (rounding & (is.null(lower) || is.null(upper)))
    stop(
      "If rounding=TRUE then you must specify the names of both the lower and upper bounds of the response.\n\n"
    )

  if (!rounding & is.null(response))
    stop("Please specify the name of the response variable.\n\n")

  if (rounding & family %in% c("negbinom","betabin"))
    stop("Rounding of responses is currently not supported for discrete response distributions. Please contact the maintainer to ask about this feature,.\n\n")

  if (is.null(jags.model.args$file))
    stop(
      "The jags.model.args list must include the name of the model file. Please see help(jags.model) for details.\n\n"
    )

  if (is.null(coda.samples.args$n.iter))
    stop(
      "The coda.samples.args list must include a variable n.iter. Please see help(coda.samples) for details.\n\n"
    )

  if(engine == "JAGS" || engine == "jags"){
    engine <- "JAGS"
    
    if (!requireNamespace("rjags", quietly = TRUE)) {
      stop("The \"rjags\" package is required to run models in JAGS. You may either install it with install.packages(\"rjags\") or run your model with \"nimble\" instead using the argument engine=\"nimble\".",
      call. = FALSE)
    }
  }
  else if(engine == "nimble"){
    if (!requireNamespace("nimble", quietly = TRUE)) {
      stop("The \"nimble\" package is required to run models in nimble. You may either install it with install.packages(\"nimble\") or run your model with \"JAGS\" instead using the argument engine=\"JAGS\".",
           call. = FALSE)
    }
    else{
      ## It is necessary to attach the namespace in order for nimble to
      ## find its own functions.
      attachNamespace("nimble")
    }
  }
  else{
    stop(engine,"is not a recognized engine for MCMC sampling.",
         call. = FALSE)
  }

  if(engine == "nimble" & family == "betabin")
    stop("The \"nimble\" package does not currently support the beta-binomial distribution. Please run your model with \"JAGS\" instead
by using the argument engine = \"JAGS\".")

  if(n.cores < 1 | !all.equal(n.cores %% 1, 0)){
    stop("Number of cores (n.cores) must be a positive integer.\n")
  }

  if(n.cores > 1){
    ## Load parallel
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("The \"parallel\" packages is required to run chains in parallel. Please install the package and try again.",
           call. = FALSE)
    }
  } 

  ## Generate JAGS input data
  message("Step 1: Generating JAGS data...")

  jags.model.args$data <-
    generateJAGSdata(
      df,
      family,
      mean.model,
      dispersion.model,
      joint.model,
      response = response,
      ntrials = ntrials,
      lower = lower,
      upper = upper,
      include.checks = include.checks,
      drop.levels = drop.levels,
      drop.missing = drop.missing
    )

  message("Done\n")

  ## Create useful variable names that will be assigned later
  mean.names.fixed <- paste0(mean.model$fixed$name,
                             ".",
                             colnames(jags.model.args$data$mean.fixed))
  
  dispersion.names.fixed <- paste0(dispersion.model$fixed$name,
                                   ".",
                                   colnames(jags.model.args$data$dispersion.fixed))

  joint.names.fixed <- paste0(joint.model$fixed$name,
                              ".",
                              colnames(jags.model.args$data$joint.fixed))

  if (!is.null(mean.model$random)) {
    mean.names.sd <-
      paste0("sd.",
             mean.model$random$name,
             ".",
             attr(terms(mean.model$random$formula), "term.labels"))

    mean.names.random <- paste0(mean.model$random$name,
                                ".",
                                colnames(jags.model.args$data$mean.random))
  }

  if (!is.null(dispersion.model$random)) {
    dispersion.names.sd <-
      paste0("sd.",
             dispersion.model$random$name,
             ".",
             attr(terms(dispersion.model$random$formula), "term.labels"))
    
    dispersion.names.random <-
      paste0(
        dispersion.model$random$name,
        ".",
        colnames(jags.model.args$data$dispersion.random)
      )
  }
  
  if (!is.null(joint.model$random)) {
      joint.names.sd <-
      paste0("sd.",
             joint.model$random$name,
             ".",
             attr(terms(joint.model$random$formula), "term.labels"))
    
    joint.names.random <-
      paste0(
        joint.model$random$name,
        ".",
        colnames(jags.model.args$data$joint.random)
      )
    }
    
  ## Perform SVD if requested
  if (svd) {
    message("    Computing singular value decompositions to improve mixing...")

    ## Compute SVD and replace original matrices with orthogonal matrices
    if(!is.null(mean.model$fixed)){
      mean.fixed.svd <- svd(jags.model.args$data$mean.fixed)
      jags.model.args$data$mean.fixed <- mean.fixed.svd$u
    }

    if(!is.null(dispersion.model$fixed)){
      dispersion.fixed.svd <-
        svd(jags.model.args$data$dispersion.fixed)
      
      jags.model.args$data$dispersion.fixed <- dispersion.fixed.svd$u
    }

    if(!is.null(joint.model$fixed)){
      joint.fixed.svd <-
        svd(jags.model.args$data$dispersion.fixed)

      jags.model.args$data$dispersion.fixed <- dispersion.fixed.svd$u
    }
        
    ## Replace design matrices with orthogal matrices

    message("Done\n")
  }

  ## Generate JAGS code
  message("Step 2: Generating JAGS code...")

  if (is.null(gencode)) {
    if (!file.exists(jags.model.args$file))
      gencode <- TRUE
    else if (overwrite)
      gencode <- TRUE
    else{
      tmp <- NULL

      while (is.null(tmp)) {
        tmp <-
          readline("The model file already exists. Do you want to overwrite it? (y/n)")

        if (tmp == "y" || tmp == "Y")
          gencode <- TRUE
        else if (tmp == "n" || tmp == "N")
          gencode <- FALSE
        else
          tmp <- NULL
      }
    }
  }

  if (gencode)
    generateJAGScode(
      family,
      jags.model.args,
      mean.model,
      dispersion.model,
      joint.model,
      rounding = rounding,
      residuals = residuals,
      include.checks = include.checks
    )

  message("Done\n")
  
  ## Generate JAGS initial values
  message("Step 3: Generating initial values...\n")

  ## Only implemented for normal model at the moment
  if(family == "gaussian"){
    if (is.null(jags.model.args$inits)) {
      if (is.null(jags.model.args$n.chains)){
        message("\n    Running three parallel chains by default...")
      }
      else {
        message("\n    Automatic generation of initial values currently works only with three chains. Setting n.chains=3...")
      }
      jags.model.args$n.chains <- 3
      
      jags.model.args$inits <-
        generateJAGSinits(family,
                          mean.model,
                          dispersion.model,
                          jags.model.args$data)
      
      message("Done\n")
    }
    else{
      if (is.null(jags.model.args$n.chains))
        jags.model.args$n.chains <-
          length(jags.model.args$inits)
      
      message("Skipped\n")
    }
  }
  else{
    jags.model.args$inits <- NULL
  }
  
  ## Save JAGS files
  if(!is.null(saveJAGSinput)){
    if(!dir.exists(saveJAGSinput))
      dir.create(saveJAGSinput)
    
    save(jags.model.args,file=file.path(saveJAGSinput,"jags_model_args.RData"))
  }

  ## Run model

  output <- list(
    df=df,
    family = family,
    mean.model = mean.model,
    dispersion.model = dispersion.model,
    joint.model = joint.model,
    jags.model.args = jags.model.args,
    coda.samples.args = coda.samples.args,
    rounding = rounding,
    parameters = parameters,
    svd = svd,
    residuals = residuals,
    drop.levels = drop.levels,
    drop.missing = drop.missing)

  if(!run.model){
    return(output)
  }

  message("Step 4: Running model\n")

  ## List parameters to monitor
  if (is.null(parameters))
    parameters <- c(
      mean.model$fixed$name,
      mean.model$random$name,
      dispersion.model$fixed$name,
      dispersion.model$random$name,
      joint.model$fixed$name,
      joint.model$random$name
    )
  
  if (residuals && !(residuals %in% parameters))
    parameters <- c(parameters, "resid")
  
  if (!is.null(mean.model$random))
    parameters <- c(parameters,
                    paste0("sd.", mean.model$random$name))
  
  if (!is.null(dispersion.model$random))
    parameters <- c(parameters,
                    paste0("sd.", dispersion.model$random$name))

      if (!is.null(joint.model$random))
        parameters <- c(parameters,
                        paste0("sd.", joint.model$random$name))
  
  coda.samples.args$variable.names <- parameters

  if(engine == "JAGS"){
    if(n.cores > 1)
      coda <- parRunJAGS(family, jags.model.args, coda.samples.args, n.cores)
    else
      coda <- runJAGS(family, jags.model.args, coda.samples.args)
  }

  if(engine == "nimble"){
    if(n.cores > 1)
      coda <- parRunNimble(jags.model.args, coda.samples.args, n.cores)
    else
      coda <- runNimble(jags.model.args, coda.samples.args)
  }

  message("Done\n")

  ## Final tidying
  message("Step 5: Tidying Output...")

  ## Identify indices of mean, dispersion, and joint parameters in coda output
  if(!is.null(mean.model$fixed))
    mean.index.fixed <- grep(paste0("^", mean.model$fixed$name),
                             coda::varnames(coda))
  else
    mean.index.fixed <- NULL

  if(!is.null(dispersion.model$fixed))
    dispersion.index.fixed <- grep(paste0("^", dispersion.model$fixed$name),
                                   coda::varnames(coda))
  else
    dispersion.index.fixed <- NULL

  if(!is.null(joint.model$fixed))
    joint.index.fixed <- grep(paste0("^", joint.model$fixed$name),
                              coda::varnames(coda))
  else
    joint.index.fixed <- NULL

  if (!is.null(mean.model$random)) {
    mean.index.sd <-
      grep(paste0("^sd\\.", mean.model$random$name), coda::varnames(coda))

    mean.index.random <-
      grep(paste0("^", mean.model$random$name),
           coda::varnames(coda))
  }

  if (!is.null(dispersion.model$random)) {
    dispersion.index.sd <-
      grep(paste0("^sd\\.", dispersion.model$random$name),
           coda::varnames(coda))

    dispersion.index.random <-
      grep(paste0("^", dispersion.model$random$name),
           coda::varnames(coda))
  }

  if (!is.null(joint.model$random)) {
    joint.index.sd <-
      grep(paste0("^sd\\.", joint.model$random$name), coda::varnames(coda))

    joint.index.random <-
      grep(paste0("^", joint.model$random$name),
           coda::varnames(coda))
  }

  ## 1) Transform chains to original scale
  if (svd) {
    for (i in 1:jags.model.args$n.chains) {
      coda[[i]][, mean.index.fixed] <-
        t(solve(mean.fixed.svd$d * t(mean.fixed.svd$v), t(coda[[i]][, mean.index.fixed])))

      coda[[i]][, dispersion.index.fixed] <-
        t(solve(
          dispersion.fixed.svd$d * t(dispersion.fixed.svd$v),
          t(coda[[i]][, dispersion.index.fixed])
        ))

      coda[[i]][, joint.index.fixed] <-
        t(solve(
          joint.fixed.svd$d * t(joint.fixed.svd$v),
          t(coda[[i]][, joint.index.fixed])
        ))
    }
  }
  
  ## Replace column names in coda with names from formula
  for (i in 1:length(coda)) {
    if(!is.null(mean.index.fixed))
       colnames(coda[[i]])[mean.index.fixed] <- mean.names.fixed

    if(!is.null(dispersion.index.fixed))
      colnames(coda[[i]])[dispersion.index.fixed] <- dispersion.names.fixed

    if(!is.null(joint.index.fixed))
      colnames(coda[[i]])[joint.index.fixed] <- joint.names.fixed

    if (!is.null(mean.model$random)) {
      colnames(coda[[i]])[mean.index.sd] <- mean.names.sd

      colnames(coda[[i]])[mean.index.random] <-
        mean.names.random
    }

    if (!is.null(dispersion.model$random)) {
      colnames(coda[[i]])[dispersion.index.sd] <- dispersion.names.sd

      colnames(coda[[i]])[dispersion.index.random] <-
        dispersion.names.random
    }

    if (!is.null(joint.model$random)) {
      colnames(coda[[i]])[joint.index.sd] <- joint.names.sd

      colnames(coda[[i]])[joint.index.random] <-
        joint.names.random
    }
  }

  message("Done\n")
  
  ## Create output object
  output$coda <- coda

  if(rounding){
    output$lower <- lower
    output$upper <- upper
  }
  else{
    output$response <- response
  }

  class(output) <- "dalmatian"

  ## Return output
  output
}

##' Prints summary information about a fitted model of class \code{dalmatian}.
##'
##' This function produces a description of the model's structure and (by default) computes and prints the summary statistics
##' computed via \code{summary.dalmatian()} and the MCMC convergence diagnostics computed via \code{convergence.dalmatian()}.
##' Further control is available by calling these functions directly.
##'
##' @title Printed Summary of a \code{dalmatian} Object
##' @param x Object of class \code{dalmatian} created by \code{dalmatian()}.
##' @param summary If TRUE (default) compute posterior summary statistics via \code{summary.dalmatian()}.
##' @param convergence If TRUE (default) compute MCMC convergence diagnostics via \code{convergence()}.
##' @param ... Ignored
##' @return List of two elements containing posterior summary statstics and convergence diagnostics (if requested). 
##' @author Simon Bonner
##' @references Bonner, S., Kim, H., Westneat, D., Mutzel, A.,
##' Wright, J., and Schofield, M.. (2021). \code{dalmatian}: A Package
##' for Fitting Double Hierarchical Linear Models in \code{R} via \code{JAGS} and
##' \code{nimble}. \emph{Journal of Statistical Software}, 100, 10, 1--25.
##' \doi{10.18637/jss.v100.i10}.
##' @export
##'
##' @examples
##' \dontrun{
##' ## Print summary of dalmatian objects
##' print(pfresults)
##' print(pfresults2)
##' }
print.dalmatian <- function(x,summary=TRUE,convergence=TRUE,...){
    ## Print information about model
    cat("Model Components\n\n")

    ## 1) Mean
    cat("  Mean:\n")
    cat("    Fixed:\n")
    cat("      Formula: ",as.character(x$mean.model$fixed$formula),"\n")
    cat("      Parameter name: ",x$mean.model$fixed$name,"\n\n")
    
    if(!is.null(x$mean.model$random)){
        cat("    Random:\n")
        cat("      Formula: ",as.character(x$mean.model$random$formula),"\n")
        cat("      Parameter name: ",x$mean.model$random$name,"\n\n")
    }

    ## 2) Dispersion
    cat("  Dispersion:\n")
    cat("    Fixed:\n")
    cat("      Formula: ",as.character(x$dispersion.model$fixed$formula),"\n")
    cat("      Parameter name: ",x$dispersion.model$fixed$name,"\n\n")
    
    if(!is.null(x$dispersion.model$random)){
        cat("    Random:\n")
        cat("      Formula: ",as.character(x$dispersion.model$random$formula),"\n")
        cat("      Parameter name: ",x$dispersion.model$random$name,"\n\n")
    }

  ## 3) Joint
  cat("  Joint:\n")
    cat("    Fixed:\n")
    cat("      Formula: ",as.character(x$joint.model$fixed$formula),"\n")
    cat("      Parameter name: ",x$joint.model$fixed$name,"\n\n")
    
    if(!is.null(x$joint.model$random)){
        cat("    Random:\n")
        cat("      Formula: ",as.character(x$joint.model$random$formula),"\n")
        cat("      Parameter name: ",x$joint.model$random$name,"\n\n")
    }

    ## Compute summary and print
    if(summary){
        summ <- summary(x)

        cat("\n\n")
        
        print(summ)
    }
    else{
        summ <- NULL
    }

    ## Compute convergence diagnostics and print
    if(convergence){
        diag <- convergence(x)

        cat("\n\nConvergence Diagnostics\n\n")

        cat("  Gelman and Rubin Diagnostics\n\n")
        print(diag$gelman)
        cat("\n\n")

        cat("  Raftery Diagnostics\n")
        print(diag$raftery)
        cat("\n\n")
    }
    else{
        diag <- NULL
    }
            
    ## Return summary and diagnostics
    list(summary=summ,
         convergence=diag)
}

##' Create traceplots and caterpillar plots from output of the fitted model.
##'
##' This function is a wrapper for the functions \code{traceplots.dalmatian()} and \code{caterpillar.dalmatian()} which
##' create traceplots and caterpillar plots of all variables stored by the sampler. Further control is available by calling
##' these functions directly.
##' 
##' @title Plot Function for \code{dalmatian} objects
##' @param x Object of class \code{dalmatian} created by \code{dalmatian()}.
##' @param trace If TRUE (default) then generate traceplots.
##' @param caterpillar If TRUE (default) then generate caterpillar plots
##' @param show If TRUE (default) then display plots as they are generated.
##' @param return_plots If TRUE (not default) return a list of \code{ggplot} objects representing the plots. 
##' @param ... Ignored
##' @return List of \code{ggplot} objects if \code{return_plots} is true. 
##' @references Bonner, S., Kim, H., Westneat, D., Mutzel, A.,
##' Wright, J., and Schofield, M.. (2021). \code{dalmatian}: A Package
##' for Fitting Double Hierarchical Linear Models in \code{R} via \code{JAGS} and
##' \code{nimble}. \emph{Journal of Statistical Software}, 100, 10, 1--25.
##' \doi{10.18637/jss.v100.i10}.
##' @export
##' @author Simon Bonner
##'
##' @examples
##' \dontrun{
##' ## Plot results for pied-flycatcher model without random effects
##' plot(pfresults)
##' 
##' ## Plot results for pied-flycatcher model with random effects
##' plot(pfresults2)
##' }
plot.dalmatian <- function(x,trace=TRUE,caterpillar=TRUE,show=TRUE,return_plots=FALSE,...){
    ## Create traceplots
    if(trace)
        traces <- traceplots(x,show=show,return_plots=return_plots)
    else
        traces <- NULL

    ## Create caterpillar plots
    if(caterpillar)
        cater <- caterpillar(x,show=show,return_plots=return_plots)
    else
        traces <- NULL

    ## Return plots if requested
    if(return_plots)
        list(traceplots=traces,caterpillar=cater)
    else
        NULL
}
