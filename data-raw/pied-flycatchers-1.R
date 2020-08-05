## ----echo=FALSE---------------------------------------------------------------
runModels <- TRUE # If true rerun models. OTW, reload previous output.


## -----------------------------------------------------------------------------
## Load package
library(dalmatian)


## -----------------------------------------------------------------------------
## Load pied flycatcher data
data(pied_flycatchers_1)


## -----------------------------------------------------------------------------
## Create variables bounding the true load
pfdata$lower=ifelse(pfdata$load==0,log(.001),log(pfdata$load-.049))
pfdata$upper=log(pfdata$load+.05)


## -----------------------------------------------------------------------------
## Mean model
mymean=list(fixed=list(name="alpha",
       formula=~ log(IVI) + broodsize + sex,
       priors=list(c("dnorm",0,.001))))

## Dispersion model
mydisp=list(fixed=list(name="psi",
      link="log",
      formula=~broodsize + sex,
      priors=list(c("dnorm",0,.001))))



## -----------------------------------------------------------------------------

## Set working directory
## By default uses a system temp directory. You probably want to change this.
workingDir <- tempdir()

## Define list of arguments for jags.model()
jm.args <- list(file=file.path(workingDir,"pied_flycatcher_1_jags.R"),n.adapt=1000)

## Define list of arguments for coda.samples()
cs.args <- list(n.iter=5000,thin=20)

## Run the model using dalmatian
## This is how the model is run. However, to save you time we will load output from a previous run instead.
if(runModels){
  pfresults <- dalmatian(df=pfdata,
                         mean.model=mymean,
                         dispersion.model=mydisp,
                         jags.model.args=jm.args,
                         coda.samples.args=cs.args,
                         rounding=TRUE,
                         lower="lower",
                         upper="upper",
                         residuals = TRUE,
                         debug=FALSE)
  
  save(pfresults,file = "pfresults.RData")
}
if(!runModels){
  load(system.file("Pied_Flycatchers_1","pfresults.RData",package="dalmatian"))
}

## -----------------------------------------------------------------------------
# Random component of mean
mymean$random=list(name="epsilon",formula=~-1 + indidx)

# Random component of dispersion
mydisp$random=list(name="xi",formula=~-1 + indidx)


## -----------------------------------------------------------------------------

## Define initial values
inits <- lapply(1:3,function(i){
                setJAGSInits(mymean,
                      mydisp,
                      y = runif(nrow(pfdata),pfdata$lower,pfdata$upper),
                      fixed.mean = tail(pfresults$coda[[i]],1)[1:4],
                      fixed.dispersion = tail(pfresults$coda[[i]],1)[5:7],
                      sd.mean = 1,
                      sd.dispersion=1)
  })

## Define list of arguments for jags.model()
jm.args <- list(file=file.path(workingDir,"pied_flycatcher_2_jags.R"),inits=inits,n.adapt=1000)

## Define list of arguments for coda.samples()
cs.args <- list(n.iter=5000,thin=10)

## Run the model using dalmatian
## This is how the model is run. However, to save you time we will load output from a previous run instead.

if(runModels){
  pfresults2 <- dalmatian(df=pfdata,
                          mean.model=mymean,
                          dispersion.model=mydisp,
                          jags.model.args=jm.args,
                          coda.samples.args=cs.args,
                          rounding=TRUE,
                          lower="lower",
                          upper="upper",
                          residuals=TRUE,
                          debug=FALSE)
  
  save(pfresults2,file = "pfresults2.RData")
}
if(!runModels){
  load(system.file("Pied_Flycatchers_1","pfresults2.RData",package="dalmatian"))
}



