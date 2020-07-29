## Load packages
library(tidyverse)

## Parameters

## Design
N <- 50 # Number of individuals
n <- 10 # Observations per individual
m <- 30  # Number of independent trials per observation

## Mean model
alpha0 <- 0
alpha1 <- 1

tau.epsilon <- 1 # Random effects variance

## Variance model
psi0 <- 0 # Intercept of variance on log scale
psi1 <- 2 # Fixed effect on dispersion

tau.xi <- 1 # Random effects variance

## Simulate data
set.seed(7777)

betabin_data_1 <- tibble(ID = as.factor(1:N),
                 x1 = rnorm(N), ## Simulate covariate for mean
                 x2 = rnorm(N), ## Simulate covariate for variance
                 z1 = rnorm(N,0,sqrt(tau.epsilon)), ## Simulate random effect for mean
                 z2 = rnorm(N,0,sqrt(tau.xi))) %>% ## Simulate random effect for variance
  crossing(tibble(Rep = 1:n,
                  m = m)) %>%
  mutate(eta.mu = alpha0 + alpha1 * x1 + z1,
         mu = (1 + exp(-eta.mu))^-1,
         eta.phi = psi0 + psi1 * x2 + z2,
         phi = (1 + exp(-eta.phi))^-1,
         alpha = mu * (1 - phi) / phi,
         beta = (1 - mu) * (1 - phi) / phi,
         p = rbeta(N * n, alpha, beta),
         y = rbinom(N * n, m, p))

## Save data to package
usethis::use_data(betabin_data_1, overwrite = TRUE)
