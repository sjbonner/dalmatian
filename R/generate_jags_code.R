generateJAGScode <- function(family,
                             jags.model.args,
                             mean.model,
                             dispersion.model,
                             joint.model,
                             rounding=FALSE,
                             residuals=FALSE){
  ## Opening header
  cat("## Created by generateJAGS: ",
      date(),"\n\n",file=jags.model.args$file)

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
  
  mean.jags <- buildLinearPredictor(component = "mean",
                                    model = mean.model$fixed,
                                    data = jags.model.args$data)
  
  ## 2b) Random effects
  if(!is.null(mean.model$random))
    mean.jags <- buildLinearPredictor(mean.jags,
                                      component = "mean",
                                      data = jags.model.args$data,
                                      model = mean.model$random,
                                      random = TRUE)

  if(!is.null(joint.model)){
    ## 2c) Joint fixed effects
    if(!is.null(joint.model$fixed))
      mean.jags <- buildLinearPredictor(mean.jags,
                                        component = "joint",
                                        data = jags.model.args$data,
                                        model = joint.model$fixed)
    
    ## 2d) Joint random effects
    if(!is.null(joint.model$random))
      mean.jags <- buildLinearPredictor(mean.jags,
                                        component = "joint",
                                        data = jags.model.args$data,
                                        model = joint.model$random,
                                        random = TRUE)
  }

  cat(mean.jags,"\n\n",file=jags.model.args$file,append=TRUE,sep="")

  ## 3) Dispersion model
  cat("\t\t ## Dispersion Model\n",file=jags.model.args$file,append=TRUE)

  ## 3a) Fixed effects
  dispersion.jags <- buildLinearPredictor(component = "dispersion",
                                          model = dispersion.model$fixed,
                                          data = jags.model.args$data)
  
  ## 3b) Random effects
  if(!is.null(dispersion.model$random))
    dispersion.jags <- buildLinearPredictor(dispersion.jags,
                                            component = "dispersion",
                                            data = jags.model.args$data,
                                            model = dispersion.model$random,
                                            random = TRUE)
  if(!is.null(joint.model)){
    ## 2c) Joint random effects
    if(!is.null(joint.model$fixed))
      dispersion.jags <- buildLinearPredictor(dispersion.jags,
                                              component = "joint",
                                              data = jags.model.args$data,
                                              model = joint.model$fixed)
    
    ## 2d) Joint fixed effects
    if(!is.null(joint.model$random))
      dispersion.jags <- buildLinearPredictor(dispersion.jags,
                                              component = "joint",
                                              data = jags.model.args$data,
                                              model = joint.model$random,
                                              random = TRUE)
  }
  
  ## Write model components to JAGS code
  cat(dispersion.jags,"\n\n",file=jags.model.args$file,append=TRUE,sep="")

  cat("\t }\n\n",file=jags.model.args$file,append=TRUE)

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

  cat("\n",file=jags.model.args$file,append=TRUE)

  if(!is.null(dispersion.model$random)){
    cat("\t ## Dispersion Model: Random\n",file=jags.model.args$file,append=TRUE)

    generatePriorsRandom(dispersion.model,jags.model.args)

    cat("\n",file=jags.model.args$file,append=TRUE)
  }

  cat("\t ## Joint Model: Fixed\n",file=jags.model.args$file,append=TRUE)
  generatePriorsFixed(joint.model,jags.model.args$file)

  if(!is.null(joint.model$random)){
    cat("\t ## Joint Model: Random\n",file=jags.model.args$file,append=TRUE)

    generatePriorsRandom(joint.model,jags.model.args)

    cat("\n",file=jags.model.args$file,append=TRUE)
  }


  ## Close model
  cat("}\n",file=jags.model.args$file,append=TRUE)
}

buildLinearPredictor <- function(lp = NULL,
                                 component,
                                 model,
                                 data,
                                 random = FALSE){

  ## If lp is empty then add lhs to linear predictor 
  if(is.null(lp)){
    if(component == "mean")
      lhs <- "muy[i]"
    else if(component == "dispersion")
      lhs <- "phi[i]"
    else
      stop("Unknown model component,", component,", in buildLinearpredictor.")

    if(!is.null(model$link))
      lp <- paste0("\t\t ",model$link,"(",lhs,") <- ")
    else 
      lp <- paste0("\t\t ", lhs," <- ")
  }
  else{
    lp <- paste0(lp,"+ ")
  }
  
  ## Extract values from model componenent
  if(random)
    p <- data[paste0(model$name,".neffects")]
  else
    p <- data[paste0(model$name,".n")]
  
  dim <- paste0("1:",p)

  ## Add predictors to lp
  matrixName <- paste0(component,".",ifelse(random,"random","fixed"))
  
  lp <- paste0(lp,
               "inprod(",
               matrixName,"[i,",dim,"],",
               model$name,"[",dim,"])")

  lp
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

