## Created by generateJAGS:  Fri Jun 12 16:36:44 2015 

model {
	 ##### Likelihood #####
	 for(i in 1:n){
		 ## Data distribution 
		 y[i] ~ dnorm(muy[i],1/pow(sdy[i],2))

		 ## Rounding
		 round1[i] <- (lower[i] < y[i])
		 round2[i] <- (y[i] < upper[i])
		 dummy[i] ~ dbern(round1[i] * round2[i])

		 ## Mean Model
		 muy[i] <- inprod(mean.fixed[i,],alpha[]) + inprod(mean.random[i,],phi1[])

		 ## Variance Model
		 log(sdy[i]) <- inprod(variance.fixed[i,],psi[])
	 }

	 ##### Priors #####
	 ## Mean Model: Fixed
	 for(k in 1:alpha.n){
		alpha[k] ~ dnorm(0,0.001)
	 }

	 ## Mean Model: Random
	 for(k in 1:phi1.n){
		phi1[k] ~ dnorm(0,tau.phi1)
	 }
	 redun.phi1~ dnorm(0,1)
	 tau.phi1~ dgamma(1.5,37.5)
	 sd.phi1<- abs(redun.phi1)/sqrt(tau.phi1)

	 ## Variance Model: Fixed
	 for(k in 1:psi.n){
		psi[k] ~ dnorm(0,0.001)
	 }
}
