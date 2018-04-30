SAS <-
function(formula, data, subset, na.action, step, model = TRUE, y = TRUE, 
	a1 = TRUE, x1 = TRUE, a2 = TRUE, x2 = TRUE, ...)
{
	call <- match.call()
	if (missing(data))
		data <- environment(formula)
	# extract the model information	
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	oformula <- as.formula(formula)
	options(warn = -1)
	formula <- as.Formula(formula)	
	if (length(formula)[2L] < 4L){
		formula <- Formula(formula(formula, rhs = 1:2))	
	}
	else {
		formula <- Formula(formula(formula, rhs = 1:4))	
	}
	mf$formula <- formula
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- terms(formula, data = data)
	mtX1 <- delete.response(terms(formula, data = data, rhs = 1L))
	mtA1 <- delete.response(terms(formula, data = data, rhs = 2L))
	mtX2 <- delete.response(terms(formula, data = data, rhs = 3L))
	mtA2 <- delete.response(terms(formula, data = data, rhs = 4L))
	# obtain response
	Y <- model.response(mf, "numeric")
	n <- length(Y)
	if (missing(step))
		step <- floor(n/log (n))
	# the first stage covariates
	X1 <- model.matrix(mtX1, mf)
	X1 <- as.matrix(X1[,-1])
	# the first stage treatment
	A1 <- model.matrix(mtA1, mf)
	A1 <- as.vector(A1[,-1])
	# the second stage covariates
	X2 <- model.matrix(mtX2, mf)
	# the second stage treatment
	A2 <- model.matrix(mtA2, mf)
	if (ncol(A2)>1){
		A2 <- as.vector(A2[,-1])
		X2 <- as.matrix(X2[,-1])
		# treatment is binary or not
		if (any((A2!=0) & (A2!=1)))
			stop("Treatment must be binary variable!")
	}
	else{
		A2 <- NULL
		X2 <- NULL
	}
	options(warn = 0)
	# treatment is binary or not
	if (any((A1!=0) & (A1!=1)))
		stop("Treatment must be binary variable!")
	if (length(Y) < 1)
		stop("empty model")
	result <- SAS.fit(y=Y, x1=X1, a1=A1, x2=X2, a2=A2, step=step)
	if (model)
		result$model <- mf
	if (y)
		result$y <- Y
	if (a1)
		result$a1 <- A1
	if (x1)
		result$x1 <- X1
	if (a2)
		result$a2 <- A2
	if (x2)
		result$x2 <- X2
	class(result) <- "SAS"
	return (result)	
}
