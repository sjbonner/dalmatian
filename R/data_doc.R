#' Pied flycatcher feeding data
#'
#' Dataset containing 5795 records of 60 pied flycatchers from 33 nest boxes feeding their nestlings during a brood manipulation experiment.
#' @format A data frame containing 5795 rows and 17 variables
#'
"pfdata"

#' Simulated data for illustrating the use of weights
#'
#' Simulated data for illustrating the use of weights in the particular case when the responses are averages of observed with different denominators
#'
#'  @format A data frame with 100 rows and 3 columns:
#'  \describe{
#'    \item{n}{The number of observations.}
#'    \item{x}{The common predictor value.}
#'    \item{y}{The mean response value.}
#'  }
#'
"weights_data_1"

#' Simulated data for illustrating the beta-binomial model
#'
#' Simulated data to show how the beta-binomial model may be fit with fixed and random effects on both the mean and dispersion.
#'
#' @format A data frame containing 500 observations and 6 columns:
#' \describe{
#'   \item{ID}{The individual ID.}
#'   \item{Rep}{The replicate number}.
#'   \item{x1}{The value of the covariate for the mean.}
#'   \item{x2}{The value of the covariate for the dispersion.}
#'   \item{m}{The number of Bernoulli trials for each observation.}
#'   \item{y}{The number of successes.}
#' }
"betabin_data_1"

' Simulated data for illustrating the negative binomial model
#'
#' Simulated data to show how the negative binomial model may be fit with fixed and random effects on both the mean and dispersion.
#'
#' @format A data frame containing 1500 observations and 5 columns:
#' \describe{
#'   \item{ID}{The individual ID.}
#'   \item{Rep}{The replicate number.}
#'   \item{x1}{The value of the covariate for the mean.}
#'   \item{x2}{The value of the covariate for the dispersion.}
#'   \item{y}{The count.}
#' }
"nbinom_data_1"

' Simulated data for illustrating the gamma model
#'
#' Simulated data to show how the gamma model may be fit with fixed and random effects on both the mean and dispersion.
#'
#' @format A data frame containing 1500 observations and 5 columns:
#' \describe{
#'   \item{ID}{The individual ID.}
#'   \item{Rep}{The replicate number.}
#'   \item{x1}{The value of the covariate for the mean.}
#'   \item{x2}{The value of the covariate for the dispersion.}
#'   \item{y}{The response.}
#' }
"gamma_data_1"
