## Reruns all of the files that create sample output for the vignettes

## Load packages
library(devtools)

## Load dalmatian
devtools::load_all()

# source("beta-binomial-1.R")

# source("negative-binomial-1.R")

# source("pied-flycatchers-1.R")

source("pied-flycatchers-2.R")

source("pied-flycatchers-joint-effects.R")
