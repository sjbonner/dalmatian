## Load packages
library(devtools)

## Load dalmatian
devtools::load_all()

# Load MCMC output
load(file.path(proj_path(),"data-mcmc","bbmcmc.RData"))

## Post-processing
bbconvergence <- convergence(bbmcmc)
bbtraceplots <- traceplots(bbmcmc, show = FALSE)
bbsummary <- summary(bbmcmc)
bbcaterpillar <- caterpillar(bbmcmc, show = FALSE)

save(bbconvergence,
     bbtraceplots,
     bbsummary,
     bbcaterpillar,
     file = file.path(proj_path(),"inst","bbresults.RData"))
