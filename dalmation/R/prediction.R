# df is a data frame

prediction <- function(object, df, method = "mean", ci = TRUE, level = 0.95) {
	
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
	## PART 3: RE-ARRANGE ESTIMATES CREATED BY DALMATION ##
	#######################################################

	# I WILL SPLIT EACH CODA CHAIN MATRIX INTO COEFFICIENT VECTORS FORM FIRST TO LAST COLUMN OF IT
	cur.index <- NULL

	# coefficients for FIXED effects in MEAN model
	mean.fixed.coef <- lapply(object$coda, function(mat) mat[,1:ncol(mean.fixed.designMat)])
	cur.index <- ncol(mean.fixed.designMat) + 1

	# coefficients for RANDOM effects in MEAN model
	if (!is.null(mean.random.designMat)) {
		mean.random.coef <- lapply(object$coda, function(mat) mat[,cur.index:(cur.index + ncol(mean.random.designMat) - 1)])
		cur.index <- cur.index + ncol(mean.random.designMat)
	}

	# coefficients for FIXED effects in VARIANCE model
	var.fixed.coef <- lapply(object$coda, function(mat) mat[,cur.index:(cur.index + ncol(var.fixed.designMat) - 1)])
	cur.index <- cur.index + ncol(var.fixed.designMat)

	# DISPERSION PARAMETER for RANDOM effects in MEAN model
	if (!is.null(mean.random.designMat)) {
		mean.disper <- lapply(object$coda, function(mat) mat[,cur.index])
		cur.index <- cur.index + 1
	}

	# DISPERSION PARAMETER and coefficients for RANDOM effects in VARIANCE model
	if (!is.null(var.random.designMat)) {
		
		var.disper <- lapply(object$coda, function(mat) mat[,cur.index])
		cur.index <- cur.index + 1

		var.random.coef <- lapply(object$coda, function(mat) mat[,cur.index:ncol(mat)])

	}

	########################################################
	## PART 4.1: PREDICTIONS FOR MEAN AND VARIANCE MODELS ##
	########################################################

	# for design matrices, duplicate them for (number of chains created) times
	mean.fixed.designList <- rep(list(mean.fixed.designMat), length(mean.fixed.coef)) # FIXED in MEAN model
	var.fixed.designList <- rep(list(var.fixed.designMat), length(var.fixed.coef)) # FIXED in VARIANCE model

	if (!is.null(mean.random.designMat)) { # RANDOM in MEAN model
		mean.random.designList <- rep(list(mean.random.designMat), length(mean.random.coef))
	}

	if(!is.null(var.random.designMat)) { # RANDOM in VARIANCE model
		var.random.designList <- rep(list(var.random.designMat), length(var.random.coef))
	}

	# FiXED effects prediction in MEAN model
	mean.fixed.pred <- Map("%*%", mean.fixed.designList, lapply(mean.fixed.coef, function(mat) t(mat)))
	# FIXED effects prediction in VARIANCE model
	var.fixed.pred <- Map("%*%", var.fixed.designList, lapply(var.fixed.coef, function(mat) t(mat)))
	
	### MEAN MODEL PREDICTION ###
	if (!is.null(mean.random.designMat)) {
		
		# RANDOM effects prediction in MEAN model
		mean.random.pred <-Map("%*%", mean.random.designList, lapply(mean.random.coef, function(mat) t(mat)))

		# SO THE FINAL PREDICTION FOR MEAN MODEL IS
		mean.pred <- Map("+", mean.fixed.pred, mean.random.pred)

	} else { mean.pred <- mean.fixed.pred }
	
	# RANDOM effects prediction in VARIANCE model
	if (!is.null(var.random.designMat)) {

		# RANDOM effects prediction in VARIANCE model
		var.random.pred <- Map("%*%", var.random.designList, lapply(var.random.coef, function(mat) t(mat)))

		# SO THE FINAL PREDICTION FOR VARIANCE MODEL is
		var.pred <- Map("+", var.fixed.pred, var.random.pred)

	} else { var.pred <- var.fixed.pred }

	
	########################################################################################
	## PART 4.2: PREDICTIONS FOR MEAN AND VARIANCE MODEL WITH MEAN (OR MODE) OF ESTIMATES ##
	########################################################################################

	if (method != "mean") { # if method == "mode"

		### get POSTERIOR MODES
		
		# mean model
		est.mean.fixed.coef <- lapply(mean.fixed.coef, function(mat) 
			apply(mat, 2, function(vec) density(vec)$x[which.max(density(vec)$y)]))
		
		if (!is.null(mean.random.designMat)) {
			
			est.mean.random.coef <- lapply(mean.random.coef, function(mat) 
				apply(mat, 2, function(vec) density(vec)$x[which.max(density(vec)$y)]))
			
			est.mean.disper <- lapply(mean.disper, function(vec) density(vec)$x[which.max(density(vec)$y)])
		}

		# variance model
		est.var.fixed.coef <- lapply(var.fixed.coef, function(mat) 
			apply(mat, 2, function(vec) density(vec)$x[which.max(density(vec)$y)]))

		if (!is.null(var.random.designMat)) {

			est.var.random.coef <- lapply(var.random.coef, function(mat) 
				apply(mat, 2, function(vec) density(vec)$x[which.max(density(vec)$y)]))

			est.var.disper <- lapply(var.disper, function(vec) density(vec)$x[which.max(density(vec)$y)])
		}

	} else {

		### get POSTERIOR MEANS
		
		# mean model
		est.mean.fixed.coef <- lapply(mean.fixed.coef, function(mat) apply(mat, 2, function(vec) mean(vec)))

		if (!is.null(mean.random.designMat)) {

			est.mean.random.coef <- lapply(mean.random.coef, function(mat) apply(mat, 2, function(vec) mean(vec)))
			est.mean.disper <- lapply(mean.disper, function(vec) mean(vec))
		}

		# variance model
		est.var.fixed.coef <- lapply(var.fixed.coef, function(mat) apply(mat, 2, function(vec) mean(vec)))

		if (!is.null(var.random.designMat)) {

			est.var.random.coef <- lapply(var.random.coef, function(mat) apply(mat, 2, function(vec) mean(vec)))
			est.var.disper <- lapply(var.disper, function(vec) mean(vec))
		}

	}

	# FiXED effects prediction in MEAN model
	est.mean.fixed.pred <- Map("%*%", mean.fixed.designList, est.mean.fixed.coef)
	# FIXED effects prediction in VARIANCE model
	est.var.fixed.pred <- Map("%*%", var.fixed.designList, est.var.fixed.coef)
	
	### MEAN MODEL PREDICTION ###
	if (!is.null(mean.random.designMat)) {
		
		# RANDOM effects prediction in MEAN model
		est.mean.random.pred <-Map("%*%", mean.random.designList, est.mean.random.coef)

		# SO THE FINAL PREDICTION FOR MEAN MODEL IS
		est.mean.pred <- Map("+", est.mean.fixed.pred, est.mean.random.pred)

	} else { est.mean.pred <- est.mean.fixed.pred }
	
	# RANDOM effects prediction in VARIANCE model
	if (!is.null(var.random.designMat)) {

		# RANDOM effects prediction in VARIANCE model
		est.var.random.pred <- Map("%*%", var.random.designList, est.var.random.coef)

		# SO THE FINAL PREDICTION FOR VARIANCE MODEL is
		est.var.pred <- Map("+", est.var.fixed.pred, est.var.random.pred)

	} else { est.var.pred <- est.var.fixed.pred }

	
	##################################################
	## PART 4.3: CREDIBLE INTERVALS FOR PREDICTIONS ##
	##################################################

	if (ci) { # if ci == TRUE

		### get CREDIBLE INTERVALS

		# mean model
		ci.mean.fixed.coef <- lapply(mean.fixed.coef, function(mat)
			apply(mat, 2, function(vec) quantile(vec, c( (1-level)/2, 1-(1 - level)/2 ))))

		if (!is.null(mean.random.designMat)) {

			ci.mean.random.coef <- lapply(mean.random.coef, function(mat)
				apply(mat, 2, function(vec) quantile(vec, c( (1-level)/2, 1-(1 - level)/2 ))))
			ci.mean.disper <- lapply(mean.disper, function(vec) quantile(vec, c( (1-level)/2, 1-(1 - level)/2 )))
		}

		# variance model
		ci.var.fixed.coef <- lapply(var.fixed.coef, function(mat)
			apply(mat, 2, function(vec) quantile(vec, c( (1-level)/2, 1-(1 - level)/2 ))))

		if (!is.null(var.random.designMat)) {

			ci.var.random.coef <- lapply(var.random.coef, function(mat)
				apply(mat, 2, function(vec) quantile(vec, c( (1-level)/2, 1-(1 - level)/2 ))))
			ci.var.disper <- lapply(var.disper, function(vec) quantile(vec, c( (1-level)/2, 1-(1 - level)/2 )))
		}

		# FiXED effects prediction in MEAN model
		ci.mean.fixed.pred <- Map("%*%", mean.fixed.designList, lapply(ci.mean.fixed.coef, function(mat) t(mat)))
		# FIXED effects prediction in VARIANCE model
		ci.var.fixed.pred <- Map("%*%", var.fixed.designList, lapply(ci.var.fixed.coef, function(mat) t(mat)))
	
		### MEAN MODEL PREDICTION ###
		if (!is.null(mean.random.designMat)) {
		
			# RANDOM effects prediction in MEAN model
			ci.mean.random.pred <-Map("%*%", mean.random.designList, lapply(ci.mean.random.coef, function(mat) t(mat)))

			# SO THE FINAL PREDICTION FOR MEAN MODEL IS
			ci.mean.pred <- Map("+", ci.mean.fixed.pred, ci.mean.random.pred)

		} else { ci.mean.pred <- ci.mean.fixed.pred }
	
		# RANDOM effects prediction in VARIANCE model
		if (!is.null(var.random.designMat)) {

			# RANDOM effects prediction in VARIANCE model
			ci.var.random.pred <- Map("%*%", var.random.designList, lapply(ci.var.random.coef, function(mat) t(mat)))

			# SO THE FINAL PREDICTION FOR VARIANCE MODEL is
			ci.var.pred <- Map("+", ci.var.fixed.pred, ci.var.random.pred)

		} else { ci.var.pred <- ci.var.fixed.pred }

	}

	
	########################################
	## PART 5: CREATE A LIST TO BE RETURN ##
	########################################

	returnList <- list()
	
	returnList$pred <- list() # predictions with posterior means or modes
	returnList$pred$mean <- est.mean.pred
	returnList$pred$var <- est.var.pred

	if (ci) {
		returnList$pred.ci <- list() # credible intervals for predictions
		returnList$pred.ci$mean <- ci.mean.pred
		returnList$pred.ci$var <- ci.var.pred
	}

	returnList$full.pred <- list() # predictions with all estimates in chains
	returnList$full.pred$mean <- mean.pred
	returnList$full.pred$var <- var.pred

	return(returnList)

}