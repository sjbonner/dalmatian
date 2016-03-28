## Load packages
library(ggplot2)
library(dhglm)
library(dalmation)
library(parallel)

## Parameter settings
pars <- list(N=100,beta=c(.02,.01),phi=2*-2.50)

## Functions
genData <- function(pars){
    ## Generate data for simple model

    ## Simulate covariate
    x <- runif(pars$N,-1,1)

    ## Compute mean and variance
    mu <- pars$beta[1] + pars$beta[2] * x
    sigmasq <- exp(pars$phi[1])

    ## Simulate response
    err <- rnorm(pars$N,rep(c(-.5,.5),rep(pars$N/2,2)))
    err <- sigmasq * err/sd(err)
    
    y <- mu + err

    ## Return data
    data.frame(x=x,y=y)
}

simFit <- function(i,pars){
    sink(file="trash")
    ## Simulate data
    data <- genData(pars)
    
    ## Fit model with DHGLM
    meanmod <- DHGLMMODELING(Model="mean",LinPred=y ~ x)
    varmod <- DHGLMMODELING(Model="dispersion",LinPred=phi ~ 1)
    
    dhglm.fit <- dhglmfit(RespDist="gaussian",DataMain=data,
                          MeanModel=meanmod,DispersionModel=varmod)

    ## Fit model with dalmation
    mymean=list(fixed=list(name="alpha",
                           formula=~ x,
                           priors=list(c("dnorm",0,.001))))
    
    myvar=list(fixed=list(name="phi",
                          link="log",
                          formula=~1,
                          priors=list(c("dnorm",0,.001))))

    jm.args <- list(file=file.path("~","Scratch","dalmation",paste0("dhglm_test_",i,".R")))

    cs.args <- list(n.adapt=100,n.iter=100,n.chain=1)

    dalmation.fit <- dalmation(df=data,
                       mean.model=mymean,
                       variance.model=myvar,
                       jags.model.args=jm.args,
                       coda.samples.args=cs.args,
                       response="y",
                       debug=FALSE)
    sink()

    list(data=data,dhglm=dhglm.fit,dalmation=dalmation.fit)
}

## Run simulations
test1 <- mclapply(1:100,simFit,mc.cores=6,pars=pars)

## Retrieve results
beta_dhglm <- t(sapply(test1,function(results) results$dhglm$beta_coeff))
phi_dhglm <- t(sapply(test1,function(results) results$dhglm$phi_coeff))
phi_dhglm2 <- phi_dhglm[,1:2]/2

summ_dalmation <- lapply(test1,function(results) summary(results$dalmation$coda))
phi_dalmation <- t(sapply(summ_dalmation,function(summ){
    c(summ[[1]][3,1:2],summ[[2]][3,c(1,5)])}
    ))

## Compute CIs
beta_dhglm <- cbind(beta_dhglm[,1:2],beta_dhglm[,1] + 1.96*outer(beta_dhglm[,2],c(-1,1)),
                  beta_dhglm[,4:5],beta_dhglm[,4] + 1.96*outer(beta_dhglm[,5],c(-1,1)))
                  
phi_dhglm <- cbind(phi_dhglm[,1:2],phi_dhglm[,1] + 1.96*outer(phi_dhglm[,2],c(-1,1)))
phi_dhglm2 <- cbind(phi_dhglm2[,1:2],phi_dhglm2[,1] + 1.96*outer(phi_dhglm2[,2],c(-1,1)))

## Compute summaries
phi_bias_dhglm <- mean(phi_dhglm[,1]) - pars$phi
phi_coverage_dhglm <- mean((phi_dhglm[,3] < pars$phi) * (phi_dhglm[,4] > pars$phi))

#phi_coverage2 <- mean((phi_dhglm2[,3] < pars$phi/2) * (phi_dhglm2[,4] > pars$phi/2))

phi_coverage_dalmation <- mean((2*phi_dalmation[,3] < pars$phi) * (2*phi_dalmation[,4] > pars$phi))

plot(phi_dhglm[,3],2*phi_dalmation[,3],pch=16)
abline(a=0,b=1)

plot(phi_dhglm[,4],2*phi_dalmation[,4],pch=16)
abline(a=0,b=1)

summary(phi_dalmation[,1:2]/(phi_dhglm[,1:2]/2))
