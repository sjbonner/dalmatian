## Load packages
library(devtools)

## Load dalmatian
devtools::load_all()

## Load MCMC output
load(file.path(proj_path(),"data-mcmc","gmcmc.RData"))

## Post-processing
gconvergence <- convergence(gmcmc)
gtraceplots <- traceplots(gmcmc, show = FALSE)
gsummary <- summary(gmcmc)
gcaterpillar <- caterpillar(gmcmc, show = FALSE)

save(gconvergence,
     gtraceplots,
     gsummary,
     gcaterpillar,
     file = file.path(proj_path(),"inst","gresults.RData"))
