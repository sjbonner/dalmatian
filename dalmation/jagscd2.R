
model {
	### likelihood ###
	for (i in 1:N) {
		y[i] ~ dnorm(mu[i], 1/pow(sigma[i], 2))

		# mean model
		mu[i] = inprod(x[i,], beta[]) + inprod(z[i,], u[])
		# variance model
		logSigma[i] = inprod(g[i,], gamma[]) + inprod(z[i,], b[])
		sigma[i] = exp(logSigma[i])
	}

	### priors ###
	# mean model: fixed effects
	for (j in 1:Nx) { beta[j] ~ dnorm(0, 1/10000) }

	# mean model: random effects
	tauMu ~ dt(0, 1, 5) I(0, )
	for (j in 1:Nz) { u[j] ~ dnorm(0, 1/pow(tauMu, 2)) }

	# variance model: fixed effects
	for (k in 1:Ng) { gamma[k] ~ dnorm(0, 1/10000)}

	# variance model: random effects
	tauSigma ~ dt(0, 1, 5) I(0, )
	for (k in 1:Nz) { b[k] ~ dnorm(0, 1/pow(tauSigma, 2)) }
}

