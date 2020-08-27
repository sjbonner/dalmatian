generateJAGScode <- function(family,
                             jags.model.args,
                             mean.model,
                             dispersion.model,
                             joint.model,
                             rounding = FALSE,
                             residuals = FALSE,
                             include.checks = TRUE){

  ## Generate header
  generateJAGSheader(jags.model.args$file)

  ## Add likelihood
  generateJAGSlhd(family,
                  jags.model.args,
                  mean.model,
                  dispersion.model,
                  joint.model,
                  rounding,
                  residuals,
                  include.checks)

  ## Add priors
  generateJAGSpriors(jags.model.args,
                     mean.model,
                     dispersion.model,
                     joint.model)

  ## Close model
  cat("}\n",file=jags.model.args$file,append=TRUE)
}

buildLinearPredictor <- function(lp = NULL,
                                 component,
                                 sep = "",
                                 model,
                                 data,
                                 random = FALSE){

  ## Add separator before next element in linear predictor
  lp <- paste0(lp,sep)
  
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
  
generateJAGSheader <- function(file){
  ## Opening header
  string <- paste("## Created by dalmatian:",date(),"\n\n")

  ## Open model
  string <- c(string,
              "model {\n")

  cat(string, file = file, sep = "")
}

generateJAGSlhd <- function(family,
                              jags.model.args,
                              mean.model,
                              dispersion.model,
                              joint.model,
                              rounding,
                              residuals,
                              include.checks = TRUE){
  
  ## Open likelihood block
  string <- paste("\t ##### Likelihood #####\n",
                  "\t for(i in 1:n){\n")

  ## Add model of mean
  string <- c(string,
              generateJAGSmeanmodel(jags.model.args,
                                    mean.model,
                                    joint.model))

  ## Add model of dispersion
  string <- c(string,
              generateJAGSdispmodel(jags.model.args,
                                    dispersion.model,
                                    joint.model))

  ## Add data model
  if(family == "gaussian")
    string <- c(string,
                generateJAGSlhd.gaussian(jags.model.args,
                                         mean.model,
                                         dispersion.model,
                                         rounding,
                                         residuals,
                                         include.checks = TRUE))

  if(family == "gamma")
    string <- c(string,
                generateJAGSlhd.gamma(jags.model.args,
                                      mean.model,
                                      dispersion.model,
                                      rounding,
                                      residuals,
                                      include.checks = TRUE))

   if(family == "betabin")
    string <- c(string,
                generateJAGSlhd.betabin(jags.model.args,
                                        mean.model,
                                        dispersion.model,
                                        rounding,
                                        residuals,
                                        include.checks = TRUE))
  
   if(family == "nbinom")
    string <- c(string,
                generateJAGSlhd.nbinom(jags.model.args,
                                       mean.model,
                                       dispersion.model,
                                       rounding,
                                       residuals,
                                       include.checks = TRUE))

  ## Close likelihood block
  string <- c(string,
              "\t }\n\n")

  ## Write string to file
  cat(string, file = jags.model.args$file, sep = "", append = TRUE)
}
  
generateJAGSrounding <- function(jags.model.args){

  ## Add code block defining rounding
  string <- c("\t\t ## Rounding\n",
              "\t\t round1[i] <- (lower[i] < y[i])\n",
              "\t\t round2[i] <- (y[i] < upper[i])\n", 
              "\t\t dummy[i] ~ dbern(round1[i] * round2[i])\n\n")

  string
}

generateJAGSchecks <- function(jags.model.args,
                               mean.check,
                               disp.check){

  ## Add comment
  string <- "\t\t ## Check range of mean and dispersion\n"
    
  ## Define probabilities based on conditions for mean and dispersion
  string <- c(string,
              "\t\t pmean.check[i] <- ", mean.check,"\n",
              "\t\t pdisp.check[i] <- ", disp.check,"\n")
  
  ## Add dummy variables with zero likelihood if conditions are false
  string <- c(string,
              "\t\t mean.check[i] ~ dbern(pmean.check[i])\n",
              "\t\t disp.check[i] ~ dbern(pdisp.check[i])\n\n")

  string
}

generateJAGSmeanmodel <- function(jags.model.args,
                                  mean.model,
                                  joint.model){
  
  ## Add comment to open block
  string <- "\t\t ## Mean Model\n"

  ## Create LHS 
  if(!is.null(mean.model$link))
    mean.lp <- paste0("\t\t ",mean.model$link,"(muy[i]) <- ")
  else 
    mean.lp <- paste0("\t\t muy[i] <- ")

  sep <- " "

  ## Create linear predictor

  ## 1) Fixed effects
  if(!is.null(mean.model$fixed)){
    mean.lp <- buildLinearPredictor(lp = mean.lp,
                                    component = "mean",
                                    sep = sep,
                                    model = mean.model$fixed,
                                    data = jags.model.args$data)

    sep <- " + "
  }
  
  ## 2) Random effects
  if(!is.null(mean.model$random)){
    mean.lp <- buildLinearPredictor(mean.lp,
                                    sep = sep,
                                    component = "mean",
                                    data = jags.model.args$data,
                                    model = mean.model$random,
                                    random = TRUE)

    sep <- " + "
  }

  ## 3) Joint effects
  if(!is.null(joint.model)){
    ## 3a) Joint fixed effects
    if(!is.null(joint.model$fixed)){
      mean.lp <- buildLinearPredictor(mean.lp,
                                      sep = sep,
                                      component = "joint",
                                      data = jags.model.args$data,
                                      model = joint.model$fixed)

      sep <- " + "
    }
    
    ## 3b) Joint random effects
    if(!is.null(joint.model$random)){
      mean.lp <- buildLinearPredictor(mean.lp,
                                      sep = sep,
                                      component = "joint",
                                      data = jags.model.args$data,
                                      model = joint.model$random,
                                      random = TRUE)

      sep <- " + "
    }
  }

  string <- c(string,
              mean.lp,
              "\n\n")

  string
}

generateJAGSdispmodel <- function(jags.model.args,
                                  dispersion.model,
                                  joint.model){
  
  ## Add comment to open code block
  string <- "\t\t ## Dispersion Model\n"

  ## Create LHS
  if(!is.null(dispersion.model$link))
    disp.lp <- paste0("\t\t ",dispersion.model$link,"(phi[i]) <- ")
  else 
    disp.lp <- paste0("\t\t phi[i] <- ")

  sep <- " "

  ## Add linear predictor
  ## 1) Fixed effects
  if(!is.null(dispersion.model$fixed)){
    disp.lp <- buildLinearPredictor(component = "dispersion",
                                    lp = disp.lp,
                                    sep = sep,
                                    model = dispersion.model$fixed,
                                    data = jags.model.args$data)

    sep <- " + "
  }
  
  ## 2) Random effects
  if(!is.null(dispersion.model$random)){
    disp.lp <- buildLinearPredictor(disp.lp,
                                    component = "dispersion",
                                    sep = sep,
                                    data = jags.model.args$data,
                                    model = dispersion.model$random,
                                    random = TRUE)

    sep <- " + "
  }

  ## 3) Joint effect
  if(!is.null(joint.model)){
    ## 3a) Joint random effects
    if(!is.null(joint.model$fixed)){
      disp.lp <- buildLinearPredictor(disp.lp,
                                      component = "joint",
                                      sep = sep,
                                      data = jags.model.args$data,
                                      model = joint.model$fixed)
      sep <- " + "
    }
  
    ## 3b) Joint fixed effects
    if(!is.null(joint.model$random)){
      disp.lp <- buildLinearPredictor(disp.lp,
                                      component = "joint",
                                      sep = sep,
                                      data = jags.model.args$data,
                                      model = joint.model$random,
                                      random = TRUE)
    }
  }

  string <- c(string,
              disp.lp,
              "\n\n")

  string
}
  
generateJAGSpriors <- function(jags.model.args,
                               mean.model,
                               dispersion.model,
                               joint.model){

  ## Open code block with comment
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

  if(!is.null(joint.model)){
    if(!is.null(joint.model$fixed)){
      cat("\t ## Joint Model: Fixed\n",file=jags.model.args$file,append=TRUE)
      generatePriorsFixed(joint.model,jags.model.args$file)
    }

    if(!is.null(joint.model$random)){
      cat("\t ## Joint Model: Random\n",file=jags.model.args$file,append=TRUE)
      
      generatePriorsRandom(joint.model,jags.model.args)
      
      cat("\n",file=jags.model.args$file,append=TRUE)
    }
  }
  
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

