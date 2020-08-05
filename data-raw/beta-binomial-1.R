## ----echo=FALSE---------------------------------------------------------------
runModels <- TRUE # If true rerun models. OTW, reload previous output.


## -----------------------------------------------------------------------------
## Load package
library(dalmatian)


## -----------------------------------------------------------------------------
## Load beta-binomial data
data(betabin_data_1)


## -----------------------------------------------------------------------------
## Define mean and variance objects
mymean <- list(fixed = list(name = "alpha",
                            formula = ~x1,
                            link = "logit",
                            priors = list(c("dnorm",0,.001))),
               random = list(name = "epsilon",
                             formula = ~ID - 1))

mydisp <- list(fixed = list(name = "psi",
                            formula = ~x2,
                            link = "logit",
                            priors = list(c("dnorm",0,.001))),
               random = list(name = "xi",
                             formula = ~ID - 1))



## -----------------------------------------------------------------------------

## Set working directory
## By default uses a system temp directory. You probably want to change this.
workingDir <- tempdir()

## Define list of arguments for jags.model()
jm.args <- list(file=file.path(workingDir,"betabin_test_1.R"),n.chains = 3, n.adapt = 1000)

## Define list of arguments for coda.samples()
cs.args <- list(n.iter=5000,thin=20)

## Run the model using dalmatian
## This is how the model is run. However, to save you time we will load output from a previous run instead.
if(runModels){
  bbresults <- dalmatian(df=betabin_data_1,
                         family = "betabin",
                         mean.model=mymean,
                         dispersion.model=mydisp,
                         jags.model.args=jm.args,
                         coda.samples.args=cs.args,
                         response = "y",
                         ntrials = "m",
                         residuals = TRUE,
                         run.model = TRUE,
                         engine = "JAGS",
                         n.cores = 3,
                         overwrite = TRUE,
                         saveJAGSinput = workingDir)
					 
  save(bbresults, file = "bbresults.RData")
}
if(!runModels){
  load(system.file("Beta_Binomial_1","bbresults.RData",package="dalmatian"))
}
