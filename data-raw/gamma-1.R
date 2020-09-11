## Load packages
library(devtools)

## Load dalmatian
devtools::load_all()

## Load negative binomial data
data(gamma_data_1)

## Define mean and variance objects
mymean <- list(fixed = list(name = "alpha",
                            formula = ~x1,
                            priors = list(c("dnorm",0,.001))),
               random = list(name = "epsilon",
                             formula = ~ID - 1),
               link = "log")

mydisp <- list(fixed = list(name = "psi",
                            formula = ~ x2,
                            priors = list(c("dnorm",0,.001))),
               random = list(name = "xi",
                             formula = ~ID - 1),
               link = "log")

## Set working directory
workingDir <- tempdir()

## Define list of arguments for jags.model()
jm.args <- list(file=file.path(workingDir,"gamma_test_1.R"),n.chains = 3, n.adapt = 1000)

## Define list of arguments for coda.samples()
cs.args <- list(n.iter=5000,thin=20)

## Run the model using dalmatian
system.time(
  gmcmc <- dalmatian(df=gamma_data_1,
                        family = "gamma",
                        mean.model=mymean,
                        dispersion.model=mydisp,
                        jags.model.args=jm.args,
                        coda.samples.args=cs.args,
                        response = "y",
                        residuals = FALSE,
                        run.model = TRUE,
                        engine = "JAGS",
                        n.cores = 3,
                        overwrite = TRUE,
                        gencode = TRUE,
                        saveJAGSinput = workingDir)
)

## For use on remote server
## save(gmcmc,file = "gmcmc.RData") 

## For use on local machine within package
save(gmcmc,
     file = file.path(proj_path(),"data-mcmc","Gamma_1","gmcmc.RData"))
