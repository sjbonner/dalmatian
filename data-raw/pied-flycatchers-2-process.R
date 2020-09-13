## Load packages
library(devtools)

## Load dalmatian
devtools::load_all()

## Load MCMC output
load(file.path(proj_path(),"data-mcmc","pfmcmc2.RData"))

## Post-processing
pfconvergence2 <- convergence(pfmcmc2)
pftraceplots2 <- traceplots(pfmcmc2, nthin = 20, show = FALSE)
pfsummary2 <- summary(pfmcmc2)
pfcaterpillar2 <- caterpillar(pfmcmc2, show = FALSE)
pfranef2 <- ranef(pfmcmc2)
pffitted2 <- predict(pfmcmc2,
                     newdata = pfdata[1:5,])

save(pfconvergence2,
     pftraceplots2,
     pfsummary2,
     pfcaterpillar2,
     pfranef2,
     pffitted2,
     file =file.path(proj_path(),"inst","pfresults2.RData"))
