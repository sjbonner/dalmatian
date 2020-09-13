## Load packages
library(devtools)

## Load dalmatian
devtools::load_all()

## Load MCMC output
load(file.path(proj_path(),"data-mcmc","pfmcmc1.RData"))

## Post-processing
pfconvergence1 <-  convergence(pfmcmc1)
pftraceplots1 <-  traceplots(pfmcmc1, show = FALSE)
pfsummary1 <-  summary(pfmcmc1)
pfcaterpillar1 <-  caterpillar(pfmcmc1, show = FALSE)
pfranef1 <-  ranef(pfmcmc1)
pffitted1 <-  predict(pfmcmc1,
                      newdata = pfdata[1:5,])

save(pfconvergence1,
     pftraceplots1,
     pfsummary1,
     pfcaterpillar1,
     pfranef1,
     pffitted1,
     file = file.path(proj_path(),"inst","pfresults1.RData"))
