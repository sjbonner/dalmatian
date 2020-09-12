## Load packages
library(devtools)

## Load dalmatian
devtools::load_all()

## Load beta-binomial data
data(betabin_data_1)

## Define mean and variance objects
mymean <- list(fixed = list(name = "alpha",
                            formula = ~x1,
                            priors = list(c("dnorm",0,.001))),
               random = list(name = "epsilon",
                             formula = ~ID - 1),
               link = "logit")

mydisp <- list(fixed = list(name = "psi",
                            formula = ~x2,
                            priors = list(c("dnorm",0,.001))),
               random = list(name = "xi",
                             formula = ~ID - 1),
               link = "logit")

## Set working directory
workingDir <- tempdir()

## Define list of arguments for jags.model()
jm.args <- list(file=file.path(workingDir,"betabin_test_1.R"),n.chains = 3, n.adapt = 1000)

## Define list of arguments for coda.samples()
cs.args <- list(n.iter=5000,thin=20)

## Run the model using dalmatian
bbmcmc <- dalmatian(df=betabin_data_1,
                       family = "betabin",
                       mean.model=mymean,
                       dispersion.model=mydisp,
                       jags.model.args=jm.args,
                       coda.samples.args=cs.args,
                       response = "y",
                       ntrials = "m",
                       residuals = FALSE,
                       run.model = TRUE,
                       engine = "JAGS",
                       n.cores = 3,
                       overwrite = TRUE,
                       saveJAGSinput = workingDir)

## For use on remote server
## save(bbmcmc,"bbmcmc.RData") 

## For use on local machine within packge
save(bbmcmc,
     file = file.path(proj_path(),"data-mcmc","bbmcmc.RData"))

## Post-processing
bbresults <- list(
  convergence = convergence(bbmcmc),
  traceplots = traceplots(bbmcmc, show = FALSE),
  summary = summary(bbmcmc),
  caterpillar = caterpillar(bbmcmc, show = FALSE))

save(bbresults,
     file = file.path(proj_path(),"inst","bbresults.RData"))
