generateJAGScode <- function(family,
                             jags.model.args,
                             mean.model,
                             dispersion.model,
                             rounding=FALSE,
                             residuals=FALSE){
  ## Opening header
  cat("## Created by generateJAGS: ", date(),"\n\n",file=jags.model.args$file)

  ## Open model
  cat("model {\n",file=jags.model.args$file,append=TRUE)

  ##### Likelihood #####
  cat("\t ##### Likelihood #####\n",file=jags.model.args$file,append=TRUE)
  cat("\t for(i in 1:n){\n",file=jags.model.args$file,append=TRUE)

  ## 1) Data model
  cat("\t\t ## Data distribution \n",file=jags.model.args$file,append=TRUE)

  if(family == "gaussian"){
    if(is.null(dispersion.model$weights))
      cat("\t\t y[i] ~ dnorm(muy[i],1/pow(phi[i],2))\n\n",file=jags.model.args$file,append=TRUE)
    else
      cat("\t\t y[i] ~ dnorm(muy[i],weights[i]/pow(phi[i],2))\n\n",file=jags.model.args$file,append=TRUE)

    cat("\t\t sdy[i] <- phi[i]\n\n",file=jags.model.args$file,append=TRUE)
  }

  else if(family == "nbinom"){
    cat("\t\t y[i] ~ dnegbin(p[i],r[i])\n",
        "\t\t r[i] <- 1/phi[i]\n",
        "\t\t p[i] <- 1/(1 + phi[i] * muy[i])\n",
        "\t\t sdy[i] <- sqrt(muy[i] / p[i])\n",
        "\n",
        file=jags.model.args$file,append=TRUE)
  }

  else if(family == "betabin"){
    cat("\t\t y[i] ~ dbetabin(alphay[i], betay[i], m[i])\n",
        "\t\t alphay[i] <- muy[i] * (1 - phi[i]) / phi[i]\n",
        "\t\t betay[i] <- (1 - muy[i]) * (1 - phi[i]) / phi[i]\n",
        "\t\t sdy[i] <- sqrt(m[i] * muy[i] * (1-muy[i]) * (1 + (m[i] - 1) * phi[i]))\n",
        "\n",
        file=jags.model.args$file,append=TRUE)
  }
  
  if(rounding){
    cat("\t\t ## Rounding\n",
        "\t\t round1[i] <- (lower[i] < y[i])\n",
        "\t\t round2[i] <- (y[i] < upper[i])\n", 
        "\t\t dummy[i] ~ dbern(round1[i] * round2[i])\n\n"
       ,file=jags.model.args$file,append=TRUE)
  }

  if(residuals){
    if(family == "betabin")
      cat("\t\t resid[i] <- (y[i] - m[i] * muy[i])/sdy[i]\n\n",
          file=jags.model.args$file,append=TRUE)
    
    else 
      cat("\t\t resid[i] <- (y[i] - muy[i])/sdy[i]\n\n",
          file=jags.model.args$file,append=TRUE)
  }

  ## 2) Mean model
  cat("\t\t ## Mean Model\n",file=jags.model.args$file,append=TRUE)
  
  ## 2a) Fixed effects
  dim.fixed <- paste0("1:",jags.model.args$data[paste0(mean.model$fixed$name,".n")])
  
  if(!is.null(mean.model$fixed$link))
    mean.jags <- paste("\t\t ",mean.model$fixed$link,"(muy[i]) <- inprod(mean.fixed[i,",dim.fixed,"],",mean.model$fixed$name,"[",dim.fixed,"])",sep="")
  else
    mean.jags <- paste("\t\t muy[i] <- inprod(mean.fixed[i,",dim.fixed,"],",mean.model$fixed$name,"[",dim.fixed,"])",sep="")
  
  ## 2b) Random effects
  if(!is.null(mean.model$random)){
    dim.random <- paste0("1:",jags.model.args$data[paste0(mean.model$random$name,".neffects")])
    
    mean.jags <- paste(mean.jags," + inprod(mean.random[i,",dim.random,"],",mean.model$random$name,"[",dim.random,"])",sep="")
  }

  cat(mean.jags,"\n\n",file=jags.model.args$file,append=TRUE,sep="")

  ## 3) Dispersion model
  cat("\t\t ## Dispersion Model\n",file=jags.model.args$file,append=TRUE)

  ## 3a) Fixed effects
   dim.fixed <- paste0("1:",jags.model.args$data[paste0(dispersion.model$fixed$name,".n")])

  if(!is.null(dispersion.model$fixed$link))
    dispersion.jags <- paste("\t\t ",dispersion.model$fixed$link,"(phi[i]) <- inprod(dispersion.fixed[i,",dim.fixed,"],",dispersion.model$fixed$name,"[",dim.fixed,"])",sep="")
  else
    dispersion.jags <- paste("\t\t phi[i] <- inprod(dispersion.fixed[i,",dim.fixed,"],",dispersion.model$fixed$name,"[",dim.fixed,"])",sep="")

  ## 3b) Random effects
  if(!is.null(dispersion.model$random)){
       dim.random <- paste0("1:",jags.model.args$data[paste0(dispersion.model$random$name,".neffects")])

       dispersion.jags <- paste(dispersion.jags," + inprod(dispersion.random[i,",dim.random,"],",dispersion.model$random$name,"[",dim.random,"])",sep="")
  }

  cat(dispersion.jags,"\n\n",file=jags.model.args$file,append=TRUE,sep="")

  ## Dispersion Components by Observation
  ## if(!is.null(mean.model$random)){
  ##     cat("\t\t ## Dispersion components by observation\n",file=jags.model.args$file,append=TRUE)

  ##     varcomp1.jags <- paste0("\t\t for(k in 1:",mean.model$random$name,".ncomponents){")
  ##     cat(varcomp1.jags,"\n",file=jags.model.args$file,append=TRUE)

  ##     varcomp2.jags <- paste0("\t\t\t varcomp.ind[i,k] <- var.",mean.model$random$name,"[k]/(sum(var.",mean.model$random$name,"[]) + sdy[i]^2)")
  ##     cat(varcomp2.jags,"\n",file=jags.model.args$file,append=TRUE)

  ##     cat(" \t\t }\n",file=jags.model.args$file,append=TRUE)

  ##     varcomp3.jags <- paste0("\t\t varcomp.ind[i,",mean.model$random$name,".ncomponents+1] <- 1 - sum(varcomp.ind[i,1:",mean.model$random$name,".ncomponents])")
  ##     cat(varcomp3.jags,"\n",file=jags.model.args$file,append=TRUE)
  ##}

  cat("\t }\n\n",file=jags.model.args$file,append=TRUE)

  ## Mean dispersion components
  ## if(!is.null(mean.model$random)){
  ##     cat("\t ## Mean dispersion components \n",file=jags.model.args$file,append=TRUE)

  ##     varcomp1.jags <- paste0("\t for(k in 1:",mean.model$random$name,".ncomponents+1){")
  ##     cat(varcomp1.jags,"\n",file=jags.model.args$file,append=TRUE)

  ##     varcomp2.jags <- paste0("\t\t varcomp[k] <- mean(varcomp.ind[,k])")
  ##     cat(varcomp2.jags,"\n",file=jags.model.args$file,append=TRUE)

  ##     cat(" \t }\n\n",file=jags.model.args$file,append=TRUE)
  ## }

  ##### Priors #####
  cat("\t ##### Priors #####\n",file=jags.model.args$file,append=TRUE)

  cat("\t ## Mean Model: Fixed\n",file=jags.model.args$file,append=TRUE)
  generatePriorsFixed(mean.model,jags.model.args$file)

  cat("\n",file=jags.model.args$file,append=TRUE)

  if(!is.null(mean.model$random)){
    cat("\t ## Mean Model: Random\n",file=jags.model.args$file,append=TRUE)
    generatePriorsRandom(mean.model,jags.model.args)

    cat("\n",file=jags.model.args$file,append=TRUE)
  }

  cat("\t ## Dispersion Model: Fixed\n",file=jags.model.args$file,append=TRUE)
  generatePriorsFixed(dispersion.model,jags.model.args$file)

  if(!is.null(dispersion.model$random)){
    cat("\t ## Dispersion Model: Random\n",file=jags.model.args$file,append=TRUE)

    generatePriorsRandom(dispersion.model,jags.model.args)

    cat("\n",file=jags.model.args$file,append=TRUE)
  }


  ## Close model
  cat("}\n",file=jags.model.args$file,append=TRUE)
}

generatePriorsFixed <- function(model,file){
  if(length(model$fixed$priors)==1){
    cat("\t for(k in 1:",model$fixed$name,".n){\n",file=file,append=TRUE,sep="")
    cat("\t\t",model$fixed$name,"[k] ~ ",model$fixed$priors[[1]][1],"(",model$fixed$priors[[1]][2],",",model$fixed$priors[[1]][3],")\n",file=file,append=TRUE,sep="")
    cat("\t }\n",file=file,append=TRUE,sep="")
  }
  else{
    for(k in 1:length(model$fixed$priors)){
      cat("\t",model$fixed$name,"[",k,"] ~ ",model$fixed$priors[[k]][1],"(",model$fixed$priors[[k]][2],",",model$fixed$priors[[k]][3],")\n",file=file,append=TRUE,sep="")
    }
  }

}

generatePriorsRandom <- function(model,jags.model.args){

  if(is.null(model$random$sigma)){
    ## Generate JAGS variable names
    redun <- paste0("redun.",model$random$name)
    redunk <- paste0(redun,"[k]")
    tau <- paste0("tau.",model$random$name,"[k]")
    sd <- paste0("sd.",model$random$name,"[k]")
    var <- paste0("var.",model$random$name,"[k]")

    ## Random effects dispersions
    cat("\t for(k in 1:",model$random$name,".ncomponents){\n",file=jags.model.args$file,append=TRUE,sep="")
    cat("\t\t ",redunk,"~ dnorm(0,1)\n",file=jags.model.args$file,append=TRUE,sep="")
    cat("\t\t ",tau,"~ dgamma(1.5,37.5)\n",file=jags.model.args$file,append=TRUE,sep="")
    #cat("\t\t ",tau,"~ dgamma(.5,.1)\n",file=jags.model.args$file,append=TRUE,sep="")
    cat("\t\t ",sd,"<- abs(",redunk,")/sqrt(",tau,")\n",file=jags.model.args$file,append=TRUE,sep="")
    cat("\t\t ",var,"<- pow(",sd,",2) \n",file=jags.model.args$file,append=TRUE,sep="")
    cat("\t }\n\n",file=jags.model.args$file,append=TRUE,sep="")

    ## Random effects
    cat("\t for(k in 1:",model$random$name,".neffects){\n",file=jags.model.args$file,append=TRUE,sep="")
    cat("\t\t",model$random$name,".tmp[k] ~ dnorm(0,tau.",model$random$name,"[",model$random$name,".levels[k]])\n",
        file=jags.model.args$file,append=TRUE,sep="")
    cat("\t\t",model$random$name,"[k] <- ",model$random$name,".tmp[k] * ",redun,"[",model$random$name,".levels[k]]\n",
        file=jags.model.args$file,append=TRUE,sep="")
    cat("\t }\n",file=jags.model.args$file,append=TRUE,sep="")
  }
  else{
    ## Random effects dispersions
    tau <- paste0("tau.",model$random$name,"[k]")
    sd <- paste0("sd.",model$random$name,"[k]")
    cat("\t for(k in 1:",model$random$name,".ncomponents){\n",file=jags.model.args$file,append=TRUE,sep="")
    cat("\t\t ",tau,"<- 1/pow(",sd,",2)","\n",file=jags.model.args$file,append=TRUE,sep="")
    cat("\t }\n\n",file=jags.model.args$file,append=TRUE,sep="")

    ## Random effects
    cat("\t for(k in 1:",model$random$name,".neffects){\n",file=jags.model.args$file,append=TRUE,sep="")
    cat("\t\t",model$random$name,"[k] ~ dnorm(0,tau.",model$random$name,"[",model$random$name,".levels[k]])\n",
        file=jags.model.args$file,append=TRUE,sep="")
    cat("\t }\n",file=jags.model.args$file,append=TRUE,sep="")
  }
}

