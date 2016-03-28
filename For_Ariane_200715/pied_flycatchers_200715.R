## Load package
library(dalmation)

## Load pied flycatcher data
data(pied_flycatchers_1)

## Create variables bounding the true load
pfdata$lower=log(pfdata$IVI-.049)
pfdata$upper=log(pfdata$IVI+.05)
