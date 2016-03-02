generateJAGScode <- function(jags.model.args,mean.model,variance.model,rounding=FALSE){
    ## Opening header
    cat("## Created by generateJAGS: ", date(),"\n\n",file=jags.model.args$file)

    ## Open model
    cat("model {\n",file=jags.model.args$file,append=TRUE)

    ##### Likelihood #####
    cat("\t ##### Likelihood #####\n",file=jags.model.args$file,append=TRUE)
    cat("\t for(i in 1:n){\n",file=jags.model.args$file,append=TRUE)

    ## 1) Data model
    cat("\t\t ## Data distribution \n",file=jags.model.args$file,append=TRUE)

    cat("\t\t y[i] ~ dnorm(muy[i],weights[i]/pow(sdy[i],2))\n\n",file=jags.model.args$file,append=TRUE)

    if(rounding){
        cat("\t\t ## Rounding\n",file=jags.model.args$file,append=TRUE)
        cat("\t\t round1[i] <- (lower[i] < y[i])\n",file=jags.model.args$file,append=TRUE)
        cat("\t\t round2[i] <- (y[i] < upper[i])\n",file=jags.model.args$file,append=TRUE)
        cat("\t\t dummy[i] ~ dbern(round1[i] * round2[i])\n\n",file=jags.model.args$file,append=TRUE)
    }

    cat("\t\t resid[i] <- (y[i] - muy[i])/sdy[i]\n\n",file=jags.model.args$file,append=TRUE)

    ## 2) Mean model
    cat("\t\t ## Mean Model\n",file=jags.model.args$file,append=TRUE)

    ## 2a) Fixed effects
    if(!is.null(mean.model$fixed$link))
        mean.jags <- paste("\t\t ",mean.model$fixed$link,"(muy[i]) <- inprod(mean.fixed[i,],",mean.model$fixed$name,"[])",sep="")
    else
        mean.jags <- paste("\t\t muy[i] <- inprod(mean.fixed[i,],",mean.model$fixed$name,"[])",sep="")

    ## 2b) Random effects
    if(!is.null(mean.model$random))
        mean.jags <- paste(mean.jags," + inprod(mean.random[i,],",mean.model$random$name,"[])",sep="")

    cat(mean.jags,"\n\n",file=jags.model.args$file,append=TRUE,sep="")

    ## 3) Variance model
    cat("\t\t ## Variance Model\n",file=jags.model.args$file,append=TRUE)

    ## 3a) Fixed effects
    if(!is.null(variance.model$fixed$link))
        variance.jags <- paste("\t\t ",variance.model$fixed$link,"(sdy[i]) <- inprod(variance.fixed[i,],",variance.model$fixed$name,"[])",sep="")
    else
        variance.jags <- paste("\t\t sdy[i] <- inprod(variance.fixed[i,],",variance.model$fixed$name,"[])",sep="")

    ## 3b) Random effects
    if(!is.null(variance.model$random))
        variance.jags <- paste(variance.jags," + inprod(variance.random[i,],",variance.model$random$name,"[])",sep="")

    cat(variance.jags,"\n\n",file=jags.model.args$file,append=TRUE,sep="")

    ## Variance Components by Observation
    ## if(!is.null(mean.model$random)){
    ##     cat("\t\t ## Variance components by observation\n",file=jags.model.args$file,append=TRUE)

    ##     varcomp1.jags <- paste0("\t\t for(k in 1:",mean.model$random$name,".ncomponents){")
    ##     cat(varcomp1.jags,"\n",file=jags.model.args$file,append=TRUE)

    ##     varcomp2.jags <- paste0("\t\t\t varcomp.ind[i,k] <- var.",mean.model$random$name,"[k]/(sum(var.",mean.model$random$name,"[]) + sdy[i]^2)")
    ##     cat(varcomp2.jags,"\n",file=jags.model.args$file,append=TRUE)

    ##     cat(" \t\t }\n",file=jags.model.args$file,append=TRUE)

    ##     varcomp3.jags <- paste0("\t\t varcomp.ind[i,",mean.model$random$name,".ncomponents+1] <- 1 - sum(varcomp.ind[i,1:",mean.model$random$name,".ncomponents])")
    ##     cat(varcomp3.jags,"\n",file=jags.model.args$file,append=TRUE)
    ##}

    cat("\t }\n\n",file=jags.model.args$file,append=TRUE)

    ## Mean variance components
    ## if(!is.null(mean.model$random)){
    ##     cat("\t ## Mean variance components \n",file=jags.model.args$file,append=TRUE)

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

    cat("\t ## Variance Model: Fixed\n",file=jags.model.args$file,append=TRUE)
    generatePriorsFixed(variance.model,jags.model.args$file)

    if(!is.null(variance.model$random)){
        cat("\t ## Variance Model: Random\n",file=jags.model.args$file,append=TRUE)

        generatePriorsRandom(variance.model,jags.model.args)

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

    ## Generate JAGS variable names
    redun <- paste0("redun.",model$random$name,"[k]")
    tau <- paste0("tau.",model$random$name,"[k]")
    sd <- paste0("sd.",model$random$name,"[k]")
    var <- paste0("var.",model$random$name,"[k]")

    ## Random effects variances
    cat("\t for(k in 1:",model$random$name,".ncomponents){\n",file=jags.model.args$file,append=TRUE,sep="")
    cat("\t\t ",redun,"~ dnorm(0,1)\n",file=jags.model.args$file,append=TRUE,sep="")
    cat("\t\t ",tau,"~ dgamma(1.5,37.5)\n",file=jags.model.args$file,append=TRUE,sep="")
    cat("\t\t ",sd,"<- abs(",redun,")/sqrt(",tau,")\n",file=jags.model.args$file,append=TRUE,sep="")
    cat("\t\t ",var,"<- pow(",sd,",2) \n",file=jags.model.args$file,append=TRUE,sep="")

    cat("\t }\n\n",file=jags.model.args$file,append=TRUE,sep="")

    ## Random effects
    cat("\t for(k in 1:",model$random$name,".neffects){\n",file=jags.model.args$file,append=TRUE,sep="")
    cat("\t\t",model$random$name,"[k] ~ dnorm(0,tau.",model$random$name,"[",model$random$name,".levels[k]])\n",
        file=jags.model.args$file,append=TRUE,sep="")
    cat("\t }\n",file=jags.model.args$file,append=TRUE,sep="")


}

