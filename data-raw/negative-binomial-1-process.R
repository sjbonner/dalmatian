## Load packages
library(devtools)

## Load dalmatian
devtools::load_all()

## Load MCMC output
load(file.path(proj_path(),"data-mcmc","nbmcmc.RData"))

## Post-processing
nbconvergence <- convergence(nbmcmc)
nbtraceplots <- traceplots(nbmcmc, show = FALSE)
nbsummary <- summary(nbmcmc)
nbcaterpillar <- caterpillar(nbmcmc, show = FALSE)

save(nbconvergence,
     nbtraceplots,
     nbsummary,
     nbcaterpillar,
     file = file.path(proj_path(),"inst","nbresults.RData"))
