##' The primary function which automates the running of \code{JAGS}.
##'
##' The primary function in the package, dalmatian automates the generation of code, data, and initial values. These are then passed as arguments to function from the \code{rjags} package which automates the generation of samplse from the posterior.
##' @title Run DGLM in \code{JAGS} via \code{rjags}
##'
##' @param df Data frame containing the response and predictor values for each individual. (data.frame)
##' @param mean.model Model list specifying the structure of the mean. (list)
##' @param variance.model Model list specifying the structure of the variance. (list)
##' @param jags.model.args  List containing named arguments of \code{jags.model}. (list)
##' @param coda.samples.args List containing named arguments of \code{coda.samples}. (list)
##' @param response Name of variable in the data frame representing the response. (character)
##' @param rounding Specifies that response has been rounded if TRUE. (logical)
##' @param lower Name of variable in the data frame representing the lower bound on the response if rounded. (character)
##' @param upper Name of variable in the data frame representing the upper bound on the response if rounded. (character)
##' @param parameters Names of parameters to monitor. If NULL then default values are selected. (character)
##' @param svd Compute Singular Variable Decomposition of model matrices to improve convergence. (logical)
##' @param debug If TRUE then enter debug model. (logical)
##' @param residuals If TRUE then compute residuals in output. (logical)
##' @param gencode If TRUE then generate code potentially overwriting existing model file. By default generate code if the file does not exist and prompt user if it does. (logical)
##' @param drop.levels If TRUE then drop unused levels from all factors in df. (logical)
##' @param drop.missing If TRUE then remove records with missing response variable. (logical)
##' @param overwrite If TRUE then overwrite existing JAGS files (non-interactive sessions only). (logical)
##' @param saveJAGSinput Directory to which jags.model input is saved prior to calling \code{jags.model()}. This is useful for debugging. No files saved if NULL. (character)
##'
##' @return samples (mcmc.list)
##' @author Simon Bonner
##' @importFrom stats terms
##' @export
dalmatian <- function(df,
                      mean.model,
                      variance.model,
                      jags.model.args,
                      coda.samples.args,
                      response = NULL,
                      rounding = FALSE,
                      lower = NULL,
                      upper = NULL,
                      parameters = NULL,
                      svd = TRUE,
                      residuals = FALSE,
                      gencode = NULL,
                      drop.levels = TRUE,
                      drop.missing = TRUE,
                      overwrite = FALSE,
                      debug = FALSE,
                      saveJAGSinput=NULL) {
  ## Enter debug state
  if (debug)
    browser()

  ## Check that input is sufficient
  if (rounding && (is.null(lower) || is.null(upper)))
    stop(
      "If rounding=TRUE then you must specify the names of both the lower and upper bounds of the response.\n\n"
    )

  if (!rounding && is.null(response))
    stop("Please specify the name of the response variable.\n\n")

  if (is.null(jags.model.args$file))
    stop(
      "The jags.model.args list must include the name of the model file. Please see help(jags.model) for details.\n\n"
    )

  if (is.null(coda.samples.args$n.iter))
    stop(
      "The coda.samples.args list must include a variable n.iter. Please see help(coda.samples) for details.\n\n"
    )

  ## Generate JAGS code
  cat("Step 1: Generating JAGS code...")

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
      jags.model.args,
      mean.model,
      variance.model,
      rounding = rounding,
      residuals = residuals
    )

  cat("Done\n")

  ## Generate JAGS input data
  cat("Step 2: Generating JAGS data...")

  jags.model.args$data <-
    generateJAGSdata(
      df,
      mean.model,
      variance.model,
      response = response,
      lower = lower,
      upper = upper,
      drop.levels = drop.levels,
      drop.missing = drop.missing
    )

  cat("Done\n")

  ## Create useful variable names that will be assigned later
  mean.names.fixed <- paste0(mean.model$fixed$name,
                             ".",
                             colnames(jags.model.args$data$mean.fixed))

  variance.names.fixed <- paste0(variance.model$fixed$name,
                                 ".",
                                 colnames(jags.model.args$data$variance.fixed))

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

  if (!is.null(variance.model$random)) {
    variance.names.sd <-
      paste0("sd.",
             variance.model$random$name,
             ".",
             attr(terms(mean.model$random$formula), "term.labels"))

    variance.names.random <-
      paste0(
        variance.model$random$name,
        ".",
        colnames(jags.model.args$data$variance.random)
      )
  }
  ## Perform SVD if requested
  if (svd) {
    cat("    Computing singular value decompositions to improve mixing...")

    ## Compute SVD
    mean.fixed.svd <- svd(jags.model.args$data$mean.fixed)
    variance.fixed.svd <-
      svd(jags.model.args$data$variance.fixed)

    ## Replace design matrices with orthogal matrices
    jags.model.args$data$mean.fixed <- mean.fixed.svd$u
    jags.model.args$data$variance.fixed <- variance.fixed.svd$u

    cat("Done\n")
  }

  ## Generate JAGS initial values
  cat("Step 3: Generating initial values...")

  if (is.null(jags.model.args$inits)) {
    if (is.null(jags.model.args$n.chains)){
      cat("\n    Running three parallel chains by default...")
    }
    else {
      cat("\n    Automatic generation of initial values currently works only with three chains. Setting n.chains=3...")
    }
    jags.model.args$n.chains <- 3

    jags.model.args$inits <-
      generateJAGSinits(mean.model,
                        variance.model,
                        jags.model.args$data)

    cat("Done\n")
  }
  else{
    if (is.null(jags.model.args$n.chains))
      jags.model.args$n.chains <-
        length(jags.model.args$inits)

    cat("Skipped\n")
  }

  ## Save JAGS files
  if(!is.null(saveJAGSinput)){
    if(!dir.exists(saveJAGSinput))
      dir.create(saveJAGSinput)
    
    save(jags.model.args,file=file.path(saveJAGSinput,"jags_model_args.RData"))
  }
  
  ## Initialize model

  cat("Step 4: Running model in JAGS\n")

  cat("    Initializing model\n")
  model <- do.call(rjags::jags.model, jags.model.args)

  ## List parameters to monitor
  if (is.null(parameters))
    parameters <- c(
      mean.model$fixed$name,
      mean.model$random$name,
      variance.model$fixed$name,
      variance.model$random$name
    )

  if (residuals && !(residuals %in% parameters))
    parameters <- c(parameters, "resid")

  if (!is.null(mean.model$random))
    parameters <- c(parameters,
                    paste0("sd.", mean.model$random$name))

  if (!is.null(variance.model$random))
    parameters <- c(parameters,
                    paste0("sd.", variance.model$random$name))

  ## Generate samples
  cat("   Generating samples\n")
  coda.samples.args$model <- model
  coda.samples.args$variable.names <- parameters

  coda <- do.call(rjags::coda.samples, coda.samples.args)
  cat("Done\n")

  ## Final tidying
  cat("Step 5: Tidying Output...")

  ## Identify indices of mean and variance parameters in coda output
  mean.index.fixed <- grep(paste0("^", mean.model$fixed$name),
                           coda::varnames(coda))

  variance.index.fixed <-
    grep(paste0("^", variance.model$fixed$name),
         coda::varnames(coda))

  if (!is.null(mean.model$random)) {
    mean.index.sd <-
      grep(paste0("^sd\\.", mean.model$random$name), coda::varnames(coda))

    mean.index.random <-
      grep(paste0("^", mean.model$random$name),
           coda::varnames(coda))
  }

  if (!is.null(variance.model$random)) {
    variance.index.sd <-
      grep(paste0("^sd\\.", variance.model$random$name), coda::varnames(coda))

    variance.index.random <-
      grep(paste0("^", variance.model$random$name),
           coda::varnames(coda))
  }


  ## 1) Transform chains to original scale
  if (svd) {
    for (i in 1:jags.model.args$n.chains) {
      coda[[i]][, mean.index.fixed] <-
        t(solve(mean.fixed.svd$d * t(mean.fixed.svd$v), t(coda[[i]][, mean.index.fixed])))

      coda[[i]][, variance.index.fixed] <-
        t(solve(
          variance.fixed.svd$d * t(variance.fixed.svd$v),
          t(coda[[i]][, variance.index.fixed])
        ))
    }
  }

  ## Replace column names in coda with names from formula
  for (i in 1:length(coda)) {
    colnames(coda[[i]])[mean.index.fixed] <- mean.names.fixed
    colnames(coda[[i]])[variance.index.fixed] <-
      variance.names.fixed

    if (!is.null(mean.model$random)) {
      colnames(coda[[i]])[mean.index.sd] <- mean.names.sd

      colnames(coda[[i]])[mean.index.random] <-
        mean.names.random
    }

    if (!is.null(variance.model$random)) {
      colnames(coda[[i]])[variance.index.sd] <- variance.names.sd

      colnames(coda[[i]])[variance.index.random] <-
        variance.names.random
    }
  }

  cat("Done\n")

  ## Create output object
  output <- list(
    mean.model = mean.model,
    variance.model = variance.model,
    jags.model.args = jags.model.args,
    coda.samples.args = coda.samples.args,
    coda = coda
  )

  class(output) <- "dalmatian"

  ## Return output
  output
}
