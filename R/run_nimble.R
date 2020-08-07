runNimble <- function(jags.model.args, coda.samples.args){

  ## Initialize model
  cat("    Initializing model\n")
  
  model <- nimble::readBUGSmodel(model = jags.model.args$file,
                                 data = jags.model.args$data,
                                 inits = jags.model.args$inits[[1]])
                           
  
  ## Generate samples
  cat("   Generating samples\n")

  coda <- nimble::nimbleMCMC(model = model,
                             inits = jags.model.args$inits,
                             monitors = coda.samples.args$variable.names,
                             niter = jags.model.args$n.adapt +
                               coda.samples.args$n.iter,
                             thin = coda.samples.args$thin,
                             nchains = jags.model.args$n.chains,
                             nburnin = jags.model.args$n.adapt,
                             samplesAsCodaMCMC = TRUE)

  return(coda)
}

parRunNimble <- function(jags.model.args, coda.samples.args, n.cores){

  ## Initialize cluster
  cl <- parallel::makeCluster(n.cores)
  
  ## Run chains
  coda.list <- parallel::parLapply(cl,
                                   X=1:jags.model.args$n.chains,
                                   fun = parRunNimble_helper,
                                   jags.model.args = jags.model.args,
                                   coda.samples.args = coda.samples.args)

  ## Combine output into a single mcmc.list object
  coda <- coda::mcmc.list(coda.list)

  return(coda)
}

parRunNimble_helper <- function(k, jags.model.args, coda.samples.args){
  ## Load nimble on node
  requireNamespace("nimble", quietly = TRUE)
  attachNamespace("nimble")
  
  ## Initialize model
  cat("    Initializing model\n")
  
  model <- nimble::readBUGSmodel(model = jags.model.args$file,
                                 data = jags.model.args$data,
                                 inits = jags.model.args$inits[[k]])
  
  ## Generate samples
  cat("   Generating samples\n")

  coda <- nimble::nimbleMCMC(model = model,
                             inits = jags.model.args$inits[[k]],
                             monitors = coda.samples.args$variable.names,
                             niter = jags.model.args$n.adapt +
                               coda.samples.args$n.iter,
                             thin = coda.samples.args$thin,
                             nchains = 1,
                             nburnin = jags.model.args$n.adapt,
                             samplesAsCodaMCMC = TRUE)

  return(coda)
}
