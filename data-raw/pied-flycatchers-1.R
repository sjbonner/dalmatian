## Load packages
library(devtools)

## Load dalmatian
devtools::load_all()

## Load pied flycatcher data
data(pied_flycatchers_1)


## Create variables bounding the true load
pfdata$lower=ifelse(pfdata$load==0,log(.001),log(pfdata$load-.049))
pfdata$upper=log(pfdata$load+.05)

##### Model 1 #####

## Mean model
mymean=list(fixed=list(name="alpha",
       formula=~ log(IVI) + broodsize + sex,
       priors=list(c("dnorm",0,.001))))

## Dispersion model
mydisp=list(fixed=list(name="psi",
      link="log",
      formula=~broodsize + sex,
      priors=list(c("dnorm",0,.001))))


## Set working directory
workingDir <- tempdir()

## Define list of arguments for jags.model()
jm.args <- list(file=file.path(workingDir,"pied_flycatcher_1_jags.R"),n.adapt=1000)

## Define list of arguments for coda.samples()
cs.args <- list(n.iter=5000,thin=20)

## Run the model using dalmatian
pfresults <- dalmatian(df=pfdata,
                       mean.model=mymean,
                       dispersion.model=mydisp,
                       jags.model.args=jm.args,
                       coda.samples.args=cs.args,
                       rounding=TRUE,
                       lower="lower",
                       upper="upper",
                       n.cores = 3,
                       residuals = FALSE,
                       overwrite = TRUE,
                       debug=FALSE)

file <- file.path(proj_path(),"inst","Pied_Flycatchers_1","pfresults.RData")
save(pfresults, file = file)

##### Model 2 #####

# Random component of mean
mymean$random=list(name="epsilon",formula=~-1 + indidx)

# Random component of dispersion
mydisp$random=list(name="xi",formula=~-1 + indidx)

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
pfresults2 <- dalmatian(df=pfdata,
                        mean.model=mymean,
                        dispersion.model=mydisp,
                        jags.model.args=jm.args,
                        coda.samples.args=cs.args,
                        rounding=TRUE,
                        lower="lower",
                        upper="upper",
                        n.cores = 3,
                        residuals=FALSE,
                        overwrite = TRUE,
                        debug=FALSE)

file <- file.path(proj_path(),"inst","Pied_Flycatchers_1","pfresults2.RData")
save(pfresults2, file = file)



