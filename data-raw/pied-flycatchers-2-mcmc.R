## Load packages
library(devtools)

## Load dalmatian
devtools::load_all()

## Load pied flycatcher data
data(pied_flycatchers_1)

## Load output from Model 1
load(file.path(proj_path(),"data-mcmc","pfmcmc1.RData"))

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
pfconvergence2 <- convergence(pfmcmc2)
pftraceplots2 <- traceplots(pfmcmc2, show = FALSE)
pfsummary2 <- summary(pfmcmc2)
pfcaterpillar2 <- caterpillar(pfmcmc2, show = FALSE)
pfranef2 <- ranef(pfmcmc2)
pffitted2 <- predict(pfmcmc2,
                     newdata = pfdata[1:5,])

save(pfconvergence2,
     pftraceplots2,
     pfsummary2,
     pfcaterpillar2,
     pfranef2,
     pffitted2,
     file =file.path(proj_path(),"inst","pfresults2.RData"))
