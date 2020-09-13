## Load packages
library(devtools)
library(tidyverse)

## Parameters

## Design
N <- 50 # Number of individuals
n <- 10 # Observations per individual
m <- 30  # Number of independent trials per observation

## Mean model
alpha0 <- 0
alpha1 <- 1

tau.epsilon <- 1 # Random effects standard deviation

## Dispersion model
psi0 <- 0 # Intercept of dispersion on log scale
psi1 <- 2 # Fixed effect on dispersion

tau.xi <- 1 # Random effects standard deviation

## Simulate data
set.seed(7777)

betabin_data_1 <- tibble(ID = as.factor(1:N),
                 x1 = rnorm(N), ## Simulate covariate for mean
                 x2 = rnorm(N), ## Simulate covariate for dipsersion
                 z1 = rnorm(N,0,tau.epsilon), ## Simulate random effect for mean
                 z2 = rnorm(N,0,tau.xi)) %>% ## Simulate random effect for dispersion
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

## Remove unobserved values
betabin_data_1 <- betabin_data_1 %>%
  select(ID, Rep, x1, x2, m, y)
         
## Save data to package
usethis::use_data(betabin_data_1, overwrite = TRUE)
