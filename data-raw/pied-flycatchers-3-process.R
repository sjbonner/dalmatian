## Load packages
library(devtools)

## Load dalmatian
devtools::load_all()

load(file.path(proj_path(),"data-mcmc","pfmcmc3.RData"))

## Post-processing
pfconvergence3 <- convergence(pfmcmc3)
pftraceplots3 <- traceplots(pfmcmc3, show = FALSE)
pfsummary3 <- summary(pfmcmc3)
pfcaterpillar3 <- caterpillar(pfmcmc3, show = FALSE)
pfranef3 <- ranef(pfmcmc3)

save(pfconvergence3,
     pftraceplots3,
     pfsummary3,
     pfcaterpillar3,
     pfranef3,
     file = file.path(proj_path(),"inst","pfresults3.RData"))
