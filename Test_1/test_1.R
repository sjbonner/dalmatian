library(rjags)

## Load PF data
data.dir <- "~/Dropbox/D_Westneat/Pied_Flycatcher/Data/"

pfdata <- read.csv(file.path(data.dir,"PFdata.csv"),header=TRUE)

## Remove observations with fecal sac removal
invalid1 <- which(pfdata$fecal==1)

## Remove boxco 5,6,7 because of odd IVI values
invalid2 <- which(pfdata$boxco %in% c(5,6,7))

## Remove observations with missing IVI or load
invalid3 <- union(which(is.na(pfdata$load)),which(is.na(pfdata$IVI)))

pfdata <- pfdata[-c(invalid1,invalid2,invalid3),]

## Renumber boxes in order
boxco <- unique(pfdata$boxco)
pfdata$boxcoidx <- sapply(pfdata$boxco,function(x) which(boxco==x))

## Renumber individuals within box in order
pfdata$indidx <- factor(pfdata$indid)

## Save data into R package
save(pfdata,file="../dalmation/data/pied_flycatchers_1.Rdata")

## Variable transformations
pfdata$logIVI <- log(pfdata$IVI)
pfdata$lower <- ifelse(pfdata$load==0,log(.001),log(pfdata$load-.049))
pfdata$upper <- log(pfdata$load+.05)

## Model definitions
mymean <- list(fixed=list(name="alpha",formula=~ logIVI + broodsize + sex,priors=list(c("dnorm",0,.001))),
               random=list(name="phi1",formula=~ indidx))

myvar <- list(fixed=list(name="psi",link="log",formula=~broodsize + sex,
                  priors=list(c("dnorm",0,.001))))

## Generate JAGS script
modelFile <- "test_1_jags.R"
generateJAGScode(modelFile,mymean,myvar,rounding=TRUE)

## Generate JAGS data
jags.data <- generateJAGSdata(pfdata,mymean,myvar,lower="lower",upper="upper")

## Generate JAGS inits
jags.inits <- generateJAGSinits(mymean,myvar,jags.data)

## Run model in JAGS
model <- jags.model(modelFile,data=jags.data,inits=jags.inits)

parameters <- c(mymean$fixed$name,mymean$random$name,myvar$fixed$name,myvar$random$name)

samples <- jags.samples(model,parameters,1000)

