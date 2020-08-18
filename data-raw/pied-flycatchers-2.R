## Load packages
library(devtools)

## Load dalmatian
devtools::load_all()

## Load pied flycatcher data
data(pied_flycatchers_1)

## Create variables bounding the true load
pfdata$lower=ifelse(pfdata$load==0,log(.001),log(pfdata$load-.049))
pfdata$upper=log(pfdata$load+.05)


# Random component of mean
mymean=list(fixed=list(name="alpha",
                formula=~ log(IVI) + broodsize + sex,
                priors=list(c("dnorm",0,.001))),
    random=list(name="epsilon",formula=~-1 + indidx + indidx:log(IVI)))

# Random component of dispersion
mydisp=list(fixed=list(name="psi",
               link="log",
               formula=~1,
               priors=list(c("dnorm",0,.001))))
## Set working directory
workingDir <- tempdir()

## Define list of arguments for jags.model()
jm.args <- list(file=file.path(workingDir,"pied_flycatcher_3_jags.R"),n.adapt=1000)

## Define list of arguments for coda.samples()
cs.args <- list(n.iter=1000)

## Run the model using dalmatian
pfresults3 <- dalmatian(df=pfdata,
                        mean.model=mymean,
                        dispersion.model=mydisp,
                        jags.model.args=jm.args,
                        coda.samples.args=cs.args,
                        rounding=TRUE,
                        lower="lower",
                        upper="upper",
                        debug=FALSE)
  
save(pfresults3,
     file = file.path(proj_path(),"Pied_Flycatchers_2","pfresults3.RData"))

