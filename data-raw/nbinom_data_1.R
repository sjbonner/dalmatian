## Load packages
library(tidyverse)

## Design
n <- 50 # Number of individuals
m <- 30  # Number of observations per individual

## Mean model
alpha0 <- log(10)
alpha1 <- log(2)

## Dispersion model
psi0 <- 0 # Intercept of dipsersion on log scale
psi1 <- 2 # Fixed effect on dispersion

tau <- 1 # Random effects variance

## Simulate data
set.seed(7777)
nbinom_data_1 <- tibble(ID = as.factor(1:n),
                 x1 = rnorm(n), ## Simulate covariate for mean
                 x2 = rnorm(n), ## Simulate covariate for dispersion
                 epsilon = rnorm(n,0,sqrt(tau)), ## Simulate random effect for mean
                 xi = rnorm(n, 0, sqrt(tau.xi))) %>% ## Simulate random effect for dispersion
mutate(eta.mu = alpha0 + alpha1 * x1 + epsilon,
       mu = exp(eta.mu),
       eta.phi = psi0 + psi1 * x2 + xi,
       phi = exp(eta.phi)) %>%
  crossing(Rep = 1:m) %>%
  mutate(y = rnbinom(n * m, 1/phi, mu = mu))

## Removed unbserved values
nbinom_data_1 <- nbinom_data_1 %>%
  select(ID, Rep, x1, x2, y)

## Save data to package
usethis::use_data(nbinom_data_1, overwrite = TRUE)
