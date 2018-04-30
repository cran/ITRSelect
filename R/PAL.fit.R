PAL.fit <-
function(y, x1, x2 = NULL, a1, a2 = NULL, IC = c("BIC", "CIC", "VIC"), 
             lambda.list = exp(seq(-3.5, 2, 0.1)), refit = TRUE, 
             control = PAL.control())
{
	n <- length(y)
	# a1 is binary or not
	if (any((a1!=0) & (a1!=1)))
		stop("Treatment must be binary variable!")
	# number of observations	
	if (length(y) < 1)
		stop("empty model")
	# concordance information criterion is used by default
	if (missing(IC))
		IC <- "CIC"
	# the penalty used in estimating the baseline and propensity score model
	penalty <- control$penalty
	if ((penalty!='SCAD')&&(penalty!='lasso')&&(penalty!='MCP'))
		stop("Undefined type of penalty function")
	
	# if it is a one-stage study
	if (is.null(a2)){
		pi2.est <- NULL
		h2.est <- NULL
		alpha2.est <- NULL
		theta2.est <- NULL
		p <- dim(as.matrix(x1))[2]
		# see if the propensity score or the baseline model is missing
		if (is.null(control$pi1.est)){
			# if missing, it is estimated via penalized logistic regression
			cat("Compute the propensity score...\n")
			alpha1 <- cv.ncvreg(x1, a1, family="binomial", penalty=penalty)
			alpha1.est <- coef(alpha1, s='lambda.min')
			alpha1.est <- alpha1.est[1:(p+1)]
			pi1.est <- plogis(as.vector(cbind(1, x1) %*% as.matrix(alpha1.est)))
			cat("Done!\n")
		}
		else{
			alpha1.est <- NULL
			pi1.est <- control$pi1.est
			if (min(pi1.est) <= 0 || max(pi1.est) >= 1)
				stop("The propensity score must be in (0,1)")
			if (length(pi1.est)==1)
				pi1.est <- rep(pi1.est, n)
		}
		if (is.null(control$h1.est)){
			# if missing, it is estimated via penalized linear regression
			cat("Compute the baseline...\n")
			theta1 <- cv.ncvreg(x1[a1==0,], y[a1==0], penalty=penalty, family="gaussian")
			theta1.est <- coef(theta1, s='lambda.min')
			theta1.est <- theta1.est[1:(p+1)]
			h1.est <- as.vector(cbind(1, x1) %*% as.matrix(theta1.est))
			cat("Done!\n")
		}
		else{
			theta1.est <- NULL
			h1.est <- control$h1.est
			if (length(h1.est)==1)
				h1.est <- rep(h1.est, n)
		}
		# apply the Dantzig selector
		if (is.null(control$kappa)){
			if (IC=="CIC"){
				kap <- 1
			}
			else if (IC=="BIC"){
				kap <- 1
			}
			else{
				kap <- 4
			}
		}
		else{
			kap <- control$kappa
		}
		cat("Compute the optimal individualized treatment regime...\n")
		beta1.est <- IC.Est.Dantzig(X=x1, Y=y, A=a1, pi.est=pi1.est, h.est=h1.est, 
									lambda.list=lambda.list, IC=IC, kap=kap, refit=refit)
		cat("Done!\n")
		beta2.est <- NULL	
	}
	# if it is a two-stage study
	else{
		# apply the Dantzig selector
		if (is.null(control$kappa)){
			if (IC=="CIC"){
				kap <- 1
			}
			else if (IC=="BIC"){
				kap <- 1
			}
			else{
				kap <- 4
			}
		}
		else{
			kap <- control$kappa
		}
		if (any((a2!=0) & (a2!=1)))
			stop("Treatment must be binary variable!")
		# see if the propensity score or baseline model is missing
		x <- cbind(x1, a1, x2)
		p <- dim(as.matrix(x))[2]
		# see if the propensity score or the baseline model is missing
		if (is.null(control$pi2.est)){
			# if missing, it is estimated via penalized logistic regression
			cat("Compute the second stage propensity score...\n")
			alpha2 <- cv.ncvreg(x, a2, family="binomial", penalty=penalty)
			alpha2.est <- coef(alpha2, s='lambda.min')
			alpha2.est <- alpha2.est[1:(p+1)]
			pi2.est <- plogis(as.vector(cbind(1, x) %*% as.matrix(alpha2.est)))
			cat("Done!\n")
		}
		else{
			alpha2.est <- NULL
			pi2.est <- control$pi2.est
			if (min(pi2.est) <= 0 || max(pi2.est) >= 1)
				stop("The propensity score must be in (0,1)")
			if (length(pi2.est)==1)
				pi2.est <- rep(pi2.est, n)
		}
		if (is.null(control$h2.est)){
			# if missing, it is estimated via penalized linear regression
			cat("Compute the second stage baseline...\n")
			theta2 <- cv.ncvreg(x[a2==0,], y[a2==0], penalty=penalty, family="gaussian")
			theta2.est <- coef(theta2, s='lambda.min')
			theta2.est <- theta2.est[1:(p+1)]
			h2.est <- as.vector(cbind(1, x) %*% as.matrix(theta2.est))
			cat("Done!\n")
		}
		else{
			theta2.est <- NULL
			h2.est <- control$h2.est
			if (length(h2.est)==1)
				h2.est <- rep(h2.est, n)
		}
		cat("Compute the second stage decision rule...\n")
		beta2.est <- IC.Est.Dantzig(X=x, Y=y, A=a2, pi.est=pi2.est, h.est=h2.est, 
									lambda.list=lambda.list, IC=IC, kap=kap, refit=refit)
		cat("Done!\n")							
		
		## estimate the pseudo outcome
		X20 <- cbind(1, x)
		if (sum(beta2.est!=0)==0)
			C2.est <- rep(0, length(y))
		else if (sum(beta2.est!=0)==1)
			C2.est <- X20[,beta2.est!=0]*beta2.est[beta2.est!=0]
		else
			C2.est <- X20[,beta2.est!=0]%*%beta2.est[beta2.est!=0]
		d2.est <- (C2.est>0)+0
		v <- y+(d2.est-a2)*C2.est
		# see if the propensity score or baseline model is missing
		p <- dim(as.matrix(x1))[2]
		# see if the propensity score or the baseline model is missing
		if (is.null(control$pi1.est)){
			# if missing, it is estimated via penalized logistic regression
			cat("Compute the first stage propensity score...\n")
			alpha1 <- cv.ncvreg(x1, a1, family="binomial", penalty=penalty)
			alpha1.est <- coef(alpha1, s='lambda.min')
			alpha1.est <- alpha1.est[1:(p+1)]
			pi1.est <- plogis(as.vector(cbind(1, x1) %*% as.matrix(alpha1.est)))
			cat("Done!\n")
		}
		else{
			alpha1.est <- NULL
			pi1.est <- control$pi1.est
			if (min(pi1.est) <= 0 || max(pi1.est) >= 1)
				stop("The propensity score must be in (0,1)")
			if (length(pi1.est)==1)
				pi1.est <- rep(pi1.est, n)
		}
		if (is.null(control$h1.est)){
			# if missing, it is estimated via penalized linear regression
			cat("Compute the first stage baseline...\n")
			theta1 <- cv.ncvreg(x1[a1==0,], v[a1==0], penalty=penalty, family="gaussian")
			theta1.est <- coef(theta1, s='lambda.min')
			theta1.est <- theta1.est[1:(p+1)]
			h1.est <- as.vector(cbind(1, x1) %*% as.matrix(theta1.est))
			cat("Done!\n")
		}
		else{
			theta1.est <- NULL
			h1.est <- control$h1.est
			if (length(h1.est)==1)
				h1.est <- rep(h1.est, n)
		}
		cat("Compute the first stage decision rule...\n")
		beta1.est <- IC.Est.Dantzig(X=x1, Y=v, A=a1, pi.est=pi1.est, h.est=h1.est, 
									lambda.list=lambda.list, IC=IC, kap=kap, refit=refit)
		cat("Done!\n")							
	}
	
	result <- list(beta2.est=beta2.est, beta1.est=beta1.est, pi2.est=pi2.est, h2.est=h2.est,
				pi1.est=pi1.est, h1.est=h1.est, alpha2.est=alpha2.est, theta2.est=theta2.est,
				alpha1.est=alpha1.est, theta1.est=theta1.est)
	class(result) <- "PAL"
	return(result)
}

## estimated concordance function
C.est <- function(beta, X, A, Y, pi.est, h.est){
  
  concor <- as.vector((A-pi.est)*(Y-h.est)/(pi.est*(1-pi.est)))
  Api <- A/pi.est
  concorApi <- as.matrix(outer(concor, Api, "*"))
  concorApi <- concorApi-t(concorApi)
  
  Xbeta <- as.vector(X%*%beta)
  IXbXb <- (outer(Xbeta, Xbeta, "-")>0)
  
  return (mean(concorApi*IXbXb))
}

## CIC criteria
CIC <- function(beta, X, A, Y, pi.est, h.est, kap=1){
  
  n <- length(Y)
  p <- dim(X)[2]
  d <- sum(abs(beta)>1e-15)
  
  return (n*C.est(beta, X, A, Y, pi.est, h.est)-d*log(p)*log(n, 10)*log(log(n, 10))/kap)
}

## VIC criteria
VIC <- function(beta, X, A, Y, pi.est, h.est, kap=4){
  
  n <- length(Y)
  p <- dim(X)[2]
  gA <- A*((cbind(1,X)%*%beta)>0)+(1-A)*((cbind(1,X)%*%beta)<=0)
  Api <- A*pi.est+(1-A)*(1-pi.est)
  gApi <- gA/Api
  d <- sum(abs(beta[2:(p+1)])>1e-15)
  
  return (n*mean(gApi*Y-(gApi-1)*(h.est+pmax(cbind(1,X)%*%beta,0)))-d*n^(1/3)*(log(p))^(2/3)*log(log(n))/kap)
}

## BIC criteria
BIC <- function(beta, X, A, Y, pi.est, h.est, kap=1){
    
  d <- sum(abs(beta)>1e-15)
  n <- length(Y)
  p <- dim(X)[2]
  X0 <- cbind(1, X)
  BIC0 <- mean(((A-pi.est)*(Y-h.est-A*(X0%*%beta)))^2)
  BIC0 <- -n*log(BIC0) - d*(log(n)+log(p+1))/kap
  
  return (BIC0)
}

## the Dantzig selector
IC.Est.Dantzig <- function(X, Y, A, pi.est, h.est, lambda.list, IC="CIC", refit=TRUE, kap){
  
  p <- dim(as.matrix(X))[2]
  lambdan <- length(lambda.list)
  
  # all the coefficients
  beta.all <- matrix(0, lambdan, p+1)
  
  # all the criteria
  obj <- rep(0, lambdan)
  
  ## standardization
  center <- colMeans(X)
  std <- sqrt(apply(X, 2, var))
  ## estimation
  X0 <- cbind(1, t((t(X)-center)/std)) 
  n <- length(A)
  
  #parameters
  flp <- c(1,rep(0,2*(p+1))) #model parameters: t,beta+,beta-
  #condition on estimating equation
  Ynew <- colSums(X0*as.vector((A-pi.est)*(Y-h.est)))
  Xnew <- crossprod(X0*as.vector(A-pi.est),(X0*as.vector(A)))
  #constraint matrix
  Alp <- Matrix(rbind(cbind(-1,-Xnew,Xnew),cbind(-1,Xnew,-Xnew),c(0,rep(1,p+1),rep(1,p+1))))
  
  cat("\n")
  for (i in 1:lambdan){    
    #constraint upper bound
    blp <- c(-Ynew,Ynew,lambda.list[i])
    
    #lower bound on variables
    dir <- rep("<=", 2*(p+1)+1)
    #lower bound on variables
    bounds <- list(lower = list(ind = 2:(2*(p+1)+1), val = rep(0, 2*(p+1))))
    result <- Rglpk_solve_LP(obj=flp, mat=Alp, dir=dir, rhs=blp, bounds=bounds, max=FALSE)
    beta.est <- result$solution[2:(p+2)]-result$solution[(p+3):(2*p+3)]
	beta.est[2:(p+1)] <- beta.est[2:(p+1)]/std
	beta.est[1] <- beta.est[1]-crossprod(center, beta.est[2:(p+1)])
    
    ###########################################
    #Double with estimating equation
    ###########################################
    if (refit&&sum(beta.est[2:(p+1)]!=0)>0&&sum(beta.est!=0)<(n/log(n))){
      index <- (1:p)[beta.est[2:(p+1)]!=0]
      X.refit <- as.matrix(cbind(1,X[, index]))
      Ynew0 <- colSums(X.refit*as.vector((A-pi.est)*(Y-h.est)))
      Xnew0 <- crossprod((X.refit*as.vector(A-pi.est)),(X.refit*as.vector(A)))
      beta.refit <- solve(Xnew0, Ynew0)
      beta.est[c(1, index+1)] <- beta.refit
    }
    
    beta.all[i, ] <- beta.est
	if (IC=="CIC")
		obj[i] <- CIC(beta.est[2:(p+1)], X, A, Y, pi.est, h.est, kap=kap)
	else if (IC=="BIC")
		obj[i] <- BIC(beta.est, X, A, Y, pi.est, h.est, kap=kap)
	else
		obj[i] <- VIC(beta.est, X, A, Y, pi.est, h.est, kap=kap)
			
	progress(100*i/lambdan)
  }
  
  return (as.vector(beta.all[which.max(obj), ]))
}
