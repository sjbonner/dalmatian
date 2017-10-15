#' Prediction Method for dalmatian Fitted Objects
#'
#' @param object Object of class \code{dalmatian} created by \code{dalmatian()}.
#' @param df data frame containing predictor values to predict response variables. (data.frame)
#' @param method Method to construct the fitted model. Either \code{"mean"} or \code{"mode"} (character)
#' @param ci returning credible intervals for predictions if TRUE (logical)
#' @param level level of credible intervals for predictions (numeric)
#' @param ... Ignored
#'
#' @return predictions (list)
#' @export
#'
predict.dalmatian <- function(object, df, method = "mean", ci = TRUE, level = 0.95,...) {

	#########################
	## PART 1: WRONG CASES ##
	#########################

	# labels for FIXED effects in mean model and variance model
	mean.fixed.label <- all.vars(object$mean.model$fixed$formula)
	var.fixed.label <- all.vars(object$variance.model$fixed$formula)

	# labels for RANDOM effects in mean model and variance model
	if (!is.null(object$mean.model$random)) { # mean model
		mean.random.label <- all.vars(object$mean.model$random$formula)
	} else {
		mean.random.label <- NULL
	}

	if (!is.null(object$variance.model$random)) { # variance model
		var.random.label <- all.vars(object$variance.model$random$formula)
	} else {
		var.random.label <- NULL
	}

	# combine all variables names
	all.label <- c(mean.fixed.label, mean.random.label, var.fixed.label, var.random.label)
	all.label <- unique(all.label)

	### CHECK IF "df" INCLUDES ALL REQUIRED VARIABLES ###
	check.names <- all.label %in% names(df)
	if (all(check.names == TRUE) == FALSE) {
        print(paste0("Missing variables: ", all.label[which(check.names == FALSE)]))
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

	### mean model predictions
	mean.fixed.pred <- mean.fixed.designMat %*% t(mean.fixed.coef)
	
	if (!is.null(mean.random.designMat)) {
	  
	  mean.random.pred <- mean.random.designMat %*% t(mean.random.coef)
	  mean.pred <- mean.fixed.pred + mean.random.pred
	  
	} else {mean.pred <- mean.fixed.pred}
	
	### variance model predictions
	var.fixed.pred <- var.fixed.designMat %*% t(var.fixed.coef)
	
	if (!is.null(var.random.designMat)) {
	  
	  var.random.pred <- var.random.designMat %*% t(var.random.coef)
	  var.pred <- var.fixed.pred + var.random.pred
	  
	} else {var.pred <- var.fixed.pred}
	
	if (method != "mean") { # if method == "mode"

		### get POSTERIOR MODES for predictions

		# mean model
		est.mean.pred <- apply(mean.pred, 1, function(vec) density(vec)$x[which.max(density(vec)$y)])

		if (!is.null(mean.random.designMat)) {
			est.mean.disper <- density(mean.disper)$x[which.max(density(mean.disper)$y)]
		}

		# variance model
		est.var.pred <- apply(var.pred, 1, function(vec) density(vec)$x[which.max(density(vec)$y)])

		if (!is.null(var.random.designMat)) {
			est.var.disper <- density(var.disper)$x[which.max(density(var.disper)$y)]
		}

	} else {

		### get POSTERIOR MEANS for predictions

		# mean model
		est.mean.pred <- apply(mean.pred, 1, function(vec) mean(vec))

		if (!is.null(mean.random.designMat)) {
			est.mean.disper <- mean(mean.disper)
		}

		# variance model
		est.var.pred <- apply(var.pred, 1, function(vec) mean(vec))

		if (!is.null(var.random.designMat)) {
			est.var.disper <- mean(var.disper)
		}

	}

	##################################################
	## PART 4.3: CREDIBLE INTERVALS FOR PREDICTIONS ##
	##################################################

	if (ci) { # if ci == TRUE

		### get CREDIBLE INTERVALS

		# mean model
		ci.mean.pred <- apply(mean.pred, 1, function(vec) quantile(vec, c( (1-level)/2, 1-(1 - level)/2 )))

		if (!is.null(mean.random.designMat)) {
			ci.mean.disper <- quantile(mean.disper, c( (1-level)/2, 1-(1 - level)/2 ))
		}

		# variance model
		ci.var.pred <- apply(var.pred, 1, function(vec) quantile(vec, c( (1-level)/2, 1-(1 - level)/2 )))

		if (!is.null(var.random.designMat)) {
			ci.var.disper <- quantile(var.disper, c( (1-level)/2, 1-(1 - level)/2 ))
		}
	}

	########################################
	## PART 5: CREATE A LIST TO BE RETURN ##
	########################################

	returnList <- list()

	mean.pred <- data.frame(Fit = est.mean.pred) # for mean model
	var.pred <- data.frame(Fit = est.var.pred) # for variance model

	if (ci) {

	  mean.pred$Lower <- ci.mean.pred[1,]
	  mean.pred$Upper <- ci.mean.pred[2,]

	  var.pred$Lower <- ci.var.pred[1,]
	  var.pred$Upper <- ci.var.pred[2,]

	}

	returnList$mean <- mean.pred
	returnList$var <- var.pred

	return(returnList)

}
