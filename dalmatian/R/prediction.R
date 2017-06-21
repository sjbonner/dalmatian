#' Prediction Method for dalmatian Fitted Objects
#'
#' @param object
#' @param df
#' @param method
#' @param ci
#' @param level
#'
#' @return
#' @export
#'
#' @examples
#'
predict.dalmatian <- function(object, df, method = "mean", ci = TRUE, level = 0.95) {

	#########################
	## PART 1: WRONG CASES ##
	#########################

	# labels for FIXED effects in mean model and variance model
	mean.fixed.label <- labels(terms(object$mean.model$fixed$formula))
	var.fixed.label <- labels(terms(object$variance.model$fixed$formula))

	# labels for RANDOM effects in mean model and variance model
	if (!is.null(object$mean.model$random)) { # mean model
		mean.random.label <- labels(terms(object$mean.model$random$formula))
	} else {
		mean.random.label <- NULL
	}

	if (!is.null(object$variance.model$random)) { # variance model
		var.random.label <- labels(terms(object$variance.model$random$formula))
	} else {
		var.random.label <- NULL
	}

	# combine all variables names
	all.label <- c(mean.fixed.label, mean.random.label, var.fixed.label, var.random.label)
	all.label <- unique(all.label)

	### CHECK IF "df" INCLUDES ALL REQUIRED VARIABLES ###
	check.names <- all.label %in% names(df)
	if (all(check.names == TRUE) == FALSE) {
		stop("df does not include all required variables. Check variable names in df.")
	}

	### CHECK IF "method" is entered correctly
	if ((method != "mean") && (method != "mode")) {
		stop("method should be either 'mean' or 'mode'.")
	}

	### CHECK IF "ci" is entered correctly
	if (!is.logical(ci)) {
		stop("ci should be a logical value: TRUE or FALSE.")
	}

	### CHECK IF "level" is entered correctly
	if (!((level > 0) && (level < 1))) {
		stop("level should be a real number between 0 and 1.")
	}

	####################################
	## PART 2: CREATE DESIGN MATIRCES ##
	####################################

	# for fixed effects in mean and variance models
	mean.fixed.designMat <- model.matrix(object$mean.model$fixed$formula, df)
	var.fixed.designMat <- model.matrix(object$variance.model$fixed$formula, df)

	# for random effects in mean and variance models
	if (!is.null(object$mean.model$random)) { # mean model
		mean.random.designMat <- model.matrix(object$mean.model$random$formula, df)
	} else {
		mean.random.designMat <- NULL
	}

	if (!is.null(object$variance.model$random)) { # variance model
		var.random.designMat <- model.matrix(object$variance.model$random$formula, df)
	} else {
		var.random.designMat <- NULL
	}

	#######################################################
	## PART 3: RE-ARRANGE ESTIMATES CREATED BY dalmatian ##
	#######################################################

	# combine all chains first
	all.chains <- do.call(rbind, object$coda)

	# I WILL SPLIT EACH CODA CHAIN MATRIX INTO COEFFICIENT VECTORS FORM FIRST TO LAST COLUMN OF IT
	cur.index <- NULL

	# coefficients for FIXED effects in MEAN model
	mean.fixed.coef <- all.chains[,1:ncol(mean.fixed.designMat)]
	cur.index <- ncol(mean.fixed.designMat) + 1

	# coefficients for RANDOM effects in MEAN model
	if (!is.null(mean.random.designMat)) {
		mean.random.coef <- all.chains[,cur.index:(cur.index + ncol(mean.random.designMat) - 1)]
		cur.index <- cur.index + ncol(mean.random.designMat)
	}

	# coefficients for FIXED effects in VARIANCE model
	var.fixed.coef <- all.chains[,cur.index:(cur.index + ncol(var.fixed.designMat) - 1)]
	cur.index <- cur.index + ncol(var.fixed.designMat)

	# DISPERSION PARAMETER for RANDOM effects in MEAN model
	if (!is.null(mean.random.designMat)) {
		mean.disper <- all.chains[,cur.index]
		cur.index <- cur.index + 1
	}

	# DISPERSION PARAMETER and coefficients for RANDOM effects in VARIANCE model
	if (!is.null(var.random.designMat)) {

		var.disper <- all.chains[,cur.index]
		cur.index <- cur.index + 1

		var.random.coef <- all.chains[,cur.index:ncol(all.chains)]

	}

	########################################################################################
	## PART 4.2: PREDICTIONS FOR MEAN AND VARIANCE MODEL WITH MEAN (OR MODE) OF ESTIMATES ##
	########################################################################################

	if (method != "mean") { # if method == "mode"

		### get POSTERIOR MODES

		# mean model
		est.mean.fixed.coef <- apply(mean.fixed.coef, 2, function(vec) density(vec)$x[which.max(density(vec)$y)])

		if (!is.null(mean.random.designMat)) {

			est.mean.random.coef <- apply(mean.random.coef, 2, function(vec) density(vec)$x[which.max(density(vec)$y)])
			est.mean.disper <- density(mean.disper)$x[which.max(density(mean.disper)$y)]
		}

		# variance model
		est.var.fixed.coef <- apply(var.fixed.coef, 2, function(vec) density(vec)$x[which.max(density(vec)$y)])

		if (!is.null(var.random.designMat)) {

			est.var.random.coef <- apply(var.random.coef, 2, function(vec) density(vec)$x[which.max(density(vec)$y)])
			est.var.disper <- density(var.disper)$x[which.max(density(var.disper)$y)]
		}

	} else {

		### get POSTERIOR MEANS

		# mean model
		est.mean.fixed.coef <- apply(mean.fixed.coef, 2, function(vec) mean(vec))

		if (!is.null(mean.random.designMat)) {

			est.mean.random.coef <- apply(mean.random.coef, 2, function(vec) mean(vec))
			est.mean.disper <- mean(mean.disper)
		}

		# variance model
		est.var.fixed.coef <- apply(var.fixed.coef, 2, function(vec) mean(vec))

		if (!is.null(var.random.designMat)) {

			est.var.random.coef <- apply(var.random.coef, 2, function(vec) mean(vec))
			est.var.disper <- mean(var.disper)
		}

	}

	# FiXED effects prediction in MEAN model
	est.mean.fixed.pred <- mean.fixed.designMat %*% est.mean.fixed.coef
	# FIXED effects prediction in VARIANCE model
	est.var.fixed.pred <- var.fixed.designMat %*% est.var.fixed.coef

	### MEAN MODEL PREDICTION ###
	if (!is.null(mean.random.designMat)) {

		# RANDOM effects prediction in MEAN model
		est.mean.random.pred <- mean.random.designMat %*% est.mean.random.coef

		# SO THE FINAL PREDICTION FOR MEAN MODEL IS
		est.mean.pred <- est.mean.fixed.pred + est.mean.random.pred

	} else { est.mean.pred <- est.mean.fixed.pred }

	# RANDOM effects prediction in VARIANCE model
	if (!is.null(var.random.designMat)) {

		# RANDOM effects prediction in VARIANCE model
		est.var.random.pred <- var.random.designMat %*% est.var.random.coef

		# SO THE FINAL PREDICTION FOR VARIANCE MODEL is
		est.var.pred <- est.var.fixed.pred + est.var.random.pred

	} else { est.var.pred <- est.var.fixed.pred }


	##################################################
	## PART 4.3: CREDIBLE INTERVALS FOR PREDICTIONS ##
	##################################################

	if (ci) { # if ci == TRUE

		### get CREDIBLE INTERVALS

		# mean model
		ci.mean.fixed.coef <- apply(mean.fixed.coef, 2, function(vec) quantile(vec, c( (1-level)/2, 1-(1 - level)/2 )))

		if (!is.null(mean.random.designMat)) {

			ci.mean.random.coef <- apply(mean.random.coef, 2, function(vec) quantile(vec, c( (1-level)/2, 1-(1 - level)/2 )))
			ci.mean.disper <- quantile(mean.disper, c( (1-level)/2, 1-(1 - level)/2 ))
		}

		# variance model
		ci.var.fixed.coef <- apply(var.fixed.coef, 2, function(vec) quantile(vec, c( (1-level)/2, 1-(1 - level)/2 )))

		if (!is.null(var.random.designMat)) {

			ci.var.random.coef <- apply(var.random.coef, 2, function(vec) quantile(vec, c( (1-level)/2, 1-(1 - level)/2 )))
			ci.var.disper <- quantile(var.disper, c( (1-level)/2, 1-(1 - level)/2 ))
		}

		# FiXED effects prediction in MEAN model
		ci.mean.fixed.pred <- mean.fixed.designMat %*% t(ci.mean.fixed.coef)
		# FIXED effects prediction in VARIANCE model
		ci.var.fixed.pred <- var.fixed.designMat %*% t(ci.var.fixed.coef)

		### MEAN MODEL PREDICTION ###
		if (!is.null(mean.random.designMat)) {

			# RANDOM effects prediction in MEAN model
			ci.mean.random.pred <- mean.random.designMat %*% t(ci.mean.random.coef)

			# SO THE FINAL PREDICTION FOR MEAN MODEL IS
			ci.mean.pred <- ci.mean.fixed.pred + ci.mean.random.pred

		} else { ci.mean.pred <- ci.mean.fixed.pred }

		# RANDOM effects prediction in VARIANCE model
		if (!is.null(var.random.designMat)) {

			# RANDOM effects prediction in VARIANCE model
			ci.var.random.pred <- var.random.designMat %*% t(ci.var.random.coef)

			# SO THE FINAL PREDICTION FOR VARIANCE MODEL is
			ci.var.pred <- ci.var.fixed.pred + ci.var.random.pred

		} else { ci.var.pred <- ci.var.fixed.pred }

	}

	########################################
	## PART 5: CREATE A LIST TO BE RETURN ##
	########################################

	returnList <- list()

	mean.pred <- data.frame(Fit = est.mean.pred) # for mean model
	var.pred <- data.frame(Fit = est.var.pred) # for variance model

	if (ci) {
	  lowerName <- paste0("Lower", level*100, "%")
	  upperName <- paste0("Upper", level*100, "%")

	  mean.pred[,2:3] <- ci.mean.pred
	  colnames(mean.pred)[2:3] <- c(lowerName, upperName)

	  var.pred[,2:3] <- ci.var.pred
	  colnames(var.pred)[2:3] <- c(lowerName, upperName)

	}

	returnList$mean <- mean.pred
	returnList$var <- var.pred

	return(returnList)

}
