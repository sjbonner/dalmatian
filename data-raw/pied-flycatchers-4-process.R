## Load packages
library(devtools)

## Load dalmatian
devtools::load_all()

load(file.path(proj_path(),"data-mcmc","pfmcmc4.RData"))

## Post-processing
pfconvergence4 <- convergence(pfmcmc4)
pftraceplots4 <- traceplots(pfmcmc4, show = FALSE)
pfsummary4 <- summary(pfmcmc4)
pfcaterpillar4 <- caterpillar(pfmcmc4, show = FALSE)
pfranef4 <- ranef(pfmcmc4)

save(pfconvergence4,
     pftraceplots4,
     pfsummary4,
     pfcaterpillar4,
     pfranef4,
     file = file.path(proj_path(),"inst","pfresults4.RData"))
