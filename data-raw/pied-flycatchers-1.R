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
                       formula=~broodsize + sex,
                       priors=list(c("dnorm",0,.001))),
            link="log")

## Set working directory
workingDir <- tempdir()

## Define list of arguments for jags.model()
jm.args <- list(file=file.path(workingDir,"pied_flycatcher_1_jags.R"),n.adapt=1000)

## Define list of arguments for coda.samples()
cs.args <- list(n.iter=5000,thin=20)

## Run the model using dalmatian
pfmcmc <- dalmatian(df=pfdata,
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

file <- file.path(proj_path(),"data-mcmc","pfmcmc.RData")
save(pfmcmc, file = file)

## Post-processing
pfresults <- list(
  convergence = convergence(pfmcmc),
  traceplots = traceplots(pfmcmc, show = FALSE),
  summary = summary(pfmcmc),
  caterpillar = caterpillar(pfmcmc, show = FALSE),
  ranef = ranef(pfmcmc),
  fitted = predict(pfmcmc,
                   newdata = pfdata[1:5,]))

save(pfresults,
     file = file.path(proj_path(),"inst","pfresults.RData"))

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
               fixed.mean = tail(pfmcmc$coda[[i]],1)[1:4],
               fixed.dispersion = tail(pfmcmc$coda[[i]],1)[5:7],
               sd.mean = 1,
               sd.dispersion=1)
})

## Define list of arguments for jags.model()
jm.args <- list(file=file.path(workingDir,"pied_flycatcher_2_jags.R"),inits=inits,n.adapt=1000)

## Define list of arguments for coda.samples()
cs.args <- list(n.iter=5000,thin=10)

## Run the model using dalmatian
pfmcmc2 <- dalmatian(df=pfdata,
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

file <- file.path(proj_path(),"data-mcmc","pfmcmc2.RData")
save(pfmcmc2, file = file)


## Post-processing
pfresults2 <- list(
  convergence = convergence(pfmcmc2),
  traceplots = traceplots(pfmcmc2, show = FALSE),
  summary = summary(pfmcmc2),
  caterpillar = caterpillar(pfmcmc2, show = FALSE),
  ranef = ranef(pfmcmc2),
  fitted = predict(pfmcmc2,
                   newdata = pfdata[1:5,]))

save(pfresults2,
     file = file.path(proj_path(),"inst","pfresults2.RData"))
