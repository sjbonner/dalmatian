## Simulate data for weighted model.

## Parameters
m = 100 # Number of groups
n = 5 + rnbinom(m,size=2,mu=10) # Observations per group

x = seq(-1,1,length=m) # Covariate

## Mean and variance
mu = 0 + 1 * x # Mean
sigmasq = exp(1 + 1 * x) # Variance

## Generate data
y = sapply(1:m,function(i) mean(rnorm(n[i],mu[i],sqrt(sigmasq[i]))))

weights.df = data.frame(n=n,x=x,y=y)

## Save data
load("weights-data-1.RData")
