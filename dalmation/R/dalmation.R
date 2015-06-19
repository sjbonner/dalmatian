##' The primary function which automates the running of \code{JAGS}.
##'
##' The primary function in the package, dalmation automates the generation of code, data, and initial values. These are then passed as arguments to function from the \code{rjags} package which automates the generation of samplse from the posterior. 
##' @title Run DGLM in \code{JAGS} via \code{rjags}
##' @param df Data frame containing the response and predictor values for each individual. (data.frame)
##' @param mean.model Model list specifying the structure of the mean. (list)
##' @param variance.model Model list specifying the structure of the variance. (list)
##' @param model.file Name of the model file. (character)
##' @param jags.model.args  List containing named arguments of \code{jags.model}. (list)
##' @param coda.samples.args List containing named arguments of \code{coda.samples}. (list)
##' @param response Name of variable in the data frame representing the response. (character)
##' @param rounding Specifies that response has been rounded if TRUE. (logical)
##' @param lower Name of variable in the data frame representing the lower bound on the response if rounded. (character)
##' @param upper Name of variable in the data frame representing the upper bound on the response if rounded. (character)
##' @param parameters Names of parameters to monitor. If NULL then default values are selected. (character)
##' @param debug If TRUE then enter debug model. (logical)
##' @return samples (mcmc.list)
##' @author Simon Bonner
##' @export
##' @importFrom coda mcmc
dalmation <- function(df,
                      mean.model,
                      variance.model,
                      jags.model.args,
                      coda.samples.args,
                      response=NULL,
                      rounding=FALSE,
                      lower=NULL,
                      upper=NULL,
                      parameters=NULL,
                      debug=FALSE){

    if(debug)
        browser()
    
    ## Check that input is sufficient
    if(rounding && (is.null(lower) || is.null(upper)))
        stop("If rounding=TRUE then you must specify the names of both the lower and upper bounds of the response.\n\n")

    if(!rounding && is.null(response))
        stop("Please specify the name of the response variable.\n\n")

    if(is.null(jags.model.args$file))
        stop("The jags.model.args list must include the name of the model file. Please see help(jags.model) for details.\n\n")

    if(is.null(coda.samples.args$n.iter))
        stop("The coda.samples.args list must include a variable n.iter. Please see help(coda.samples) for details.\n\n")
    
    ## Generate JAGS code
    cat("Step 1: Generating JAGS code...")
    
    generateJAGScode(jags.model.args$file,mymean,myvar,rounding=rounding)

    cat("Done\n")

    ## Generate JAGS input data
    cat("Step 2: Generating JAGS data...")
    
    jags.model.args$data <- generateJAGSdata(df,mean.model,variance.model,lower=lower,upper=upper)

    cat("Done\n")

    ## Generate JAGS initial values
    cat("Step 3: Generating initial values...")

    if(is.null(jags.model.args$inits)){
        if(is.null(jags.model.args$n.chains))
            cat("\n    Running three parallel chains by default...")
        jags.model.args$n.chains <- 3
        
        jags.model.args$inits <- generateJAGSinits(mean.model,variance.model,jags.model.args$data,jags.model.args$n.chains)
        
        cat("Done\n")
    }
    else{
        if(is.null(jags.model.args$n.chains))
            jags.model.args$n.chains <- length(jags.model.args$inits)
            
        cat("Skipped\n")
    }

    ## Initialize model

    cat("Step 4: Running model in JAGS\n")

    cat("    Initializing model\n")
    model <- do.call(rjags::jags.model,jags.model.args)
    
    ## List parameters to monitor
    if(is.null(parameters))
        parameters <- c(mean.model$fixed$name,
                        mean.model$random$name,
                        variance.model$fixed$name,
                        variance.model$random$name)

    if(!is.null(mean.model$random))
        parameters <- c(parameters,
                        paste0("sd.",mean.model$random$name))

    if(!is.null(variance.model$random))
        parameters <- c(parameters,
                        paste0("sd.",variance.model$random$name))
    
    ## Generate 1000 samples
    cat("   Generating samples\n")
    coda.samples.args$model <- model
    coda.samples.args$variable.names <- parameters
    
    coda <- do.call(rjags::coda.samples,coda.samples.args)
    cat("Done\n")
    
    ## Return list of output
    return(list(jags.model.args=jags.model.args,
                coda.samples.args=coda.samples.args,
                coda=coda))
}
