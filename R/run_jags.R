runJAGS <- function(family,jags.model.args, coda.samples.args){

  ## Load modules
  rjags::load.module("glm")
  
  if(family == "betabin")
    rjags::load.module("mix")
  
  ## Initialize model
  cat("    Initializing model\n")
  model <- do.call(rjags::jags.model, jags.model.args)

  ## Generate samples
  cat("   Generating samples\n")
  coda.samples.args$model <- model
  coda <- do.call(rjags::coda.samples, coda.samples.args)

  return(coda)
}

parRunJAGS <- function(family,jags.model.args, coda.samples.args, n.cores){

  ## Load dclone
  if (!requireNamespace("dclone", quietly = TRUE)) {
    stop("The \"dclone\" packages is required to run chains in parallel. Please install the package and try again.",
         call. = FALSE)
  }

  ## Initialize cluster
  cl <- parallel::makeCluster(n.cores)
  
  jags.model.args$cl <-
    coda.samples.args$cl <- cl
    
  jags.model.args$name <- "dalmatian"

  ## Load modules
  parallel::clusterCall(cl,rjags::load.module,"glm")
  
  if(family == "betabin")
    parallel::clusterCall(cl,rjags::load.module,"mix")
  
  ## Initialize model
  cat("    Initializing model\n")
  do.call(dclone::parJagsModel, jags.model.args)

  ## Generate samples
  cat("   Generating samples\n")
  coda.samples.args$cl <- cl
  coda.samples.args$model <- "dalmatian"
  coda <- do.call(dclone::parCodaSamples, coda.samples.args)
    
  ## Close cluster
  parallel::stopCluster(cl)

  return(coda)
}
  
