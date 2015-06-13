generateJAGScode <- function(file,mean.model,variance.model,priors,rounding=FALSE){
    ## Opening header
    cat("## Created by generateJAGS: ", date(),"\n\n",file=file)

    ## Open model
    cat("model {\n",file=file,append=TRUE)

    ##### Likelihood #####
    cat("\t ##### Likelihood #####\n",file=file,append=TRUE)
    cat("\t for(i in 1:n){\n",file=file,append=TRUE)
        
    ## 1) Data model
    cat("\t\t ## Data distribution \n",file=file,append=TRUE)
    cat("\t\t y[i] ~ dnorm(muy[i],1/pow(sdy[i],2))\n\n",file=file,append=TRUE)

    if(rounding){
        cat("\t\t ## Rounding\n",file=file,append=TRUE)
        cat("\t\t round1[i] <- (lower[i] < y[i])\n",file=file,append=TRUE)
        cat("\t\t round2[i] <- (y[i] < upper[i])\n",file=file,append=TRUE)
        cat("\t\t dummy[i] ~ dbern(round1[i] * round2[i])\n\n",file=file,append=TRUE)
    }

    ## 2) Mean model
    cat("\t\t ## Mean Model\n",file=file,append=TRUE)

    if(!is.null(mean.model$fixed$link))
        mean.jags <- paste("\t\t ",mean.model$fixed$link,"(muy[i]) <- inprod(mean.fixed[i,],",mean.model$fixed$name,"[])",sep="")
    else
        mean.jags <- paste("\t\t muy[i] <- inprod(mean.fixed[i,],",mean.model$fixed$name,"[])",sep="")

    if(!is.null(mean.model$random))
        mean.jags <- paste(mean.jags," + inprod(mean.random[i,],",mean.model$random$name,"[])",sep="")

    cat(mean.jags,"\n\n",file=file,append=TRUE,sep="")
    
    ## 3) Variance model
    cat("\t\t ## Variance Model\n",file=file,append=TRUE)

    if(!is.null(variance.model$fixed$link))
        variance.jags <- paste("\t\t ",variance.model$fixed$link,"(sdy[i]) <- inprod(variance.fixed[i,],",variance.model$fixed$name,"[])\n",sep="")
    else
        variance.jags <- paste("\t\t sdy[i] <- inprod(variance.fixed[i,],",variance.model$fixed$name,"[])\n\n",sep="")

    cat(variance.jags,file=file,append=TRUE,sep="")

    cat("\t }\n\n",file=file,append=TRUE)
    
    ##### Priors #####
    cat("\t ##### Priors #####\n",file=file,append=TRUE)

    cat("\t ## Mean Model: Fixed\n",file=file,append=TRUE)
    generatePriorsFixed(mean.model,file)

    cat("\n",file=file,append=TRUE)

    if(!is.null(mean.model$random)){
        cat("\t ## Mean Model: Random\n",file=file,append=TRUE)
        generatePriorsRandom(mean.model,file)
        
        cat("\n",file=file,append=TRUE)
    }
        
    cat("\t ## Variance Model: Fixed\n",file=file,append=TRUE)
    generatePriorsFixed(variance.model,file)

    if(!is.null(variance.model$random)){
        cat("\t ## Variance Model: Random\n",file=file,append=TRUE)
        
        generatePriorsRandom(variance.model,file)
        
        cat("\n",file=file,append=TRUE)
    }

    ## Close model
    cat("}\n",file=file,append=TRUE)
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

generatePriorsRandom <- function(model,file){

    redun <- paste0("redun.",model$random$name)
    tau <- paste0("tau.",model$random$name)
    sd <- paste0("sd.",model$random$name)
    
    cat("\t for(k in 1:",model$random$name,".n){\n",file=file,append=TRUE,sep="")
    
    cat("\t\t",model$random$name,"[k] ~ dnorm(0,tau.",model$random$name,")\n",file=file,append=TRUE,sep="")
    cat("\t }\n",file=file,append=TRUE,sep="")
    cat("\t ",redun,"~ dnorm(0,1)\n",file=file,append=TRUE,sep="")
    cat("\t ",tau,"~ dgamma(1.5,37.5)\n",file=file,append=TRUE,sep="")
    cat("\t ",sd,"<- abs(",redun,")/sqrt(",tau,")\n",file=file,append=TRUE,sep="")
}

