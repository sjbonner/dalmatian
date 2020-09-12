## Load packages
library(devtools)

## Load dalmatian
devtools::load_all()

## Load pied flycatcher data
data(pied_flycatchers_1)

## Create variables bounding the true load
pfdata$lower=ifelse(pfdata$load==0,log(.001),log(pfdata$load-.049))
pfdata$upper=log(pfdata$load+.05)

## Mean model
mymean=list(fixed=list(name="alpha",
                       formula=~ log(IVI) + sex,
                       priors=list(c("dnorm",0,.001))))

## Dispersion model
mydisp=list(fixed=list(name="psi",
                       formula=~sex,
                       priors=list(c("dnorm",0,.001))),
            link="log")

## Joint components
myjoint = list(fixed = list(name = "gamma",
                            formula = ~-1 + broodsize,
                            priors = list(c("dnorm",0,.001))))

## Set working directory
## By default uses a system temp directory. You probably want to change this.
workingDir <- tempdir()

## Define list of arguments for jags.model()
jm.args <- list(file=file.path(workingDir,"pied_flycatcher_joint_1_jags.R"),n.adapt=1000)

## Define list of arguments for coda.samples()
cs.args <- list(n.iter=1000,thin=20)

## Run the model using dalmatian
pfresults4 <- dalmatian(df=pfdata,
                        mean.model=mymean,
                        dispersion.model=mydisp,
                        joint.model = myjoint,
                        jags.model.args=jm.args,
                        coda.samples.args=cs.args,
                        rounding=TRUE,
                        lower="lower",
                        upper="upper",
                        residuals = FALSE,
                        run.model = TRUE,
                        n.cores = 3,
                        overwrite = TRUE,
                        debug=FALSE)
  
file <- file.path(proj_path(),"inst","Pied_Flycatchers_3","pfresults4.RData")
save(pfresults4, file = file)

