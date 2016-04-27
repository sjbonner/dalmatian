##' Upgrade dalmation object
##'
##' Upgrades dalmation object for use with new postprocessing functions.
##'
##' @param input Existing object to upgrade
##'
##' @title Upgrade dalmation object
##'
##' @return output (dalmation)
##' @author Simon Bonner
##' @export
upgradeDalmation <- function(input) {

  ## Attach input
  attach(input)

  ## Create useful variable names that will be assigned later
  if (!is.null(mean.model$random)) {
    mean.names.sd <-
      paste0("sd.",
             mean.model$random$name,
             ".",
             attr(terms(mean.model$random$formula), "term.labels"))
 }

  if (!is.null(variance.model$random)) {
    variance.names.sd <-
      paste0("sd.",
             variance.model$random$name,
             ".",
             attr(terms(mean.model$random$formula), "term.labels"))
 }

  ## Identify indices of mean and variance parameters in coda output
  if (!is.null(mean.model$random)) {
    mean.index.sd <-
      grep(paste0("^sd\\.", mean.model$random$name), coda::varnames(coda))
  }

  if (!is.null(variance.model$random)) {
    variance.index.sd <-
      grep(paste0("^sd\\.", variance.model$random$name), coda::varnames(coda))
  }

  ## Replace column names in coda with names from formula
  for (i in 1:length(coda)) {
     if (!is.null(mean.model$random)) {
      colnames(coda[[i]])[mean.index.sd] <- mean.names.sd
    }

    if (!is.null(variance.model$random)) {
      colnames(coda[[i]])[variance.index.sd] <- variance.names.sd
    }
  }

  cat("Done\n")

  ## Create output object
  output <- list(
    mean.model = mean.model,
    variance.model = variance.model,
    jags.model.args = jags.model.args,
    coda.samples.args = coda.samples.args,
    coda = coda
  )

  ## Change class
  class(output) <- "dalmation"

  ## Detach input
  detach(input)

  ## Return output
  output
}
