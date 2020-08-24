## Load packages
library(tidyverse)

## Design
n <- 30 # Number of individuals
m <- 10  # Number of observations per individual

## Mean model
alpha0 <- log(10)
alpha1 <- log(2)

tau.epsilon <- 1 # Random effects standard deviation

## Dispersion model
psi0 <- 0 # Intercept of dipsersion on log scale
psi1 <- 0 #2 # Fixed effect on dispersion

tau.xi <- 1 # Random effects standard deviation

## Simulate data
set.seed(7777)
gamma_data_1 <- tibble(ID = as.factor(1:n),
                 x1 = rnorm(n), ## Simulate covariate for mean
                 x2 = rnorm(n), ## Simulate covariate for dispersion
                 epsilon = rnorm(n,0,sqrt(tau.epsilon)), ## Simulate random effect for mean
                 xi = 0) %>% #rnorm(n, 0, sqrt(tau.xi))) %>% ## Simulate random effect for dispersion
mutate(eta.mu = alpha0 + alpha1 * x1 + epsilon,
       mu = exp(eta.mu),
       eta.phi = psi0 + psi1 * x2 + xi,
       phi = exp(eta.phi)) %>%
  crossing(Rep = 1:m) %>%
  mutate(y = rgamma(n * m, shape = mu/phi, rate = 1/phi))

## Removed unbserved values
gamma_data_1 <- gamma_data_1 %>%
  select(ID, Rep, x1, x2, epsilon, xi, y)

## Save data to package
usethis::use_data(gamma_data_1, overwrite = TRUE)
