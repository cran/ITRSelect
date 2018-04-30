SAS.fit <-
function(y, x1, x2 = NULL, a1, a2 = NULL, step)
{
	n <- length(y)
	if (missing(step))
		step <- floor(n/log(n))
	# a1 is binary or not
	if (any((a1!=0) & (a1!=1)))
		stop("Treatment must be binary variable!")
	# number of observations	
	if (length(y) < 1)
		stop("empty model")
	
	# if it is a one-stage study
	if (is.null(a2)){
		# see if the propensity score or the baseline model is missing
		cat("Compute the optimal individualized treatment regime...\n")
		beta1.est <- ISS(X=x1, Y=y, A=a1, step=step)$beta
		cat("Done!\n")
		beta2.est <- NULL	
	}
	# if it is a two-stage study
	else{
		# see if treatments are binary
		if (any((a2!=0) & (a2!=1)))
			stop("Treatment must be binary variable!")
		# see if the propensity score or baseline model is missing
		x <- cbind(x1, a1, x2)
		cat("Compute the second stage decision rule...\n")
		beta2.est <- ISS(X=x, Y=y, A=a2, step=step)
		v <- beta2.est$V
		beta2.est <- beta2.est$beta
		cat("Done!\n")							
		
		## estimate the pseudo outcome
		cat("Compute the first stage decision rule...\n")
		beta1.est <- ISS(X=x1, Y=v, A=a1, step=step)$beta
		cat("Done!\n")							
	}
	
	result <- list(beta2.est=beta2.est, beta1.est=beta1.est)
	class(result) <- "SAS"
	return(result)
}


ISS <- function(X, Y, A, step){
  #some constants
  n <- length(Y)
  p <- dim(X)[2]
  
  # V_est: use the Q-learning method to estimate
  V_est<-matrix(0, nrow=n, ncol=step+1)
  a_opt <- matrix(0, nrow=n, ncol=(step+1))
  S <- rep(0,step+1)
  regime_coef <- matrix(0, nrow=(step+1), ncol=(p+1))
  
  reg0 <- lm(Y~A)
  a_opt[,1] <- as.numeric(reg0$coefficients[2]>=0)
  regime_coef[1,] <- c(as.numeric(reg0$coefficients[2]>=0),rep(0,p))
  V_est[,1] <- reg0$coef[1]+reg0$coef[2]*a_opt[,1]
  S[1] <- sum(V_est[,1]-Y)
  
  #select the most important variable
  #S score
  S_temp <- rep(0,p)
  coef_all <- matrix(0,nrow=p,ncol=2)
  V_temp <- matrix(0,nrow=n,ncol=p)
  cat("\n")
  for(j in 1:p){
    reg <- lm(Y~X[,j]*A)
    coef <- reg$coef
    coef_all[j,] <- coef[3:4]
    diff <- coef[3]+coef[4]*X[,j]
    g_opt <- as.numeric(diff>=0)
    S_temp[j] <- sum(diff*(g_opt-a_opt[,1]))
	V_temp[,j] <- Y+diff*(g_opt-A)
  }
  progress(100/step)
  #the most important variable
  imp <- which.max(S_temp)
  coef_temp <- coef_all[imp,]
  regime_coef[2,c(1,imp+1)] <- coef_temp
  a_opt[,2] <- as.numeric(cbind(1,X[,imp])%*%coef_temp>0)
  V_est[,2] <- V_temp[,imp]
  S[2] <- S_temp[imp]
  
  #select the rest important variables sequentially
  for (j in 2:step){
    S_temp <- rep(0,p) #in the jth search
    coef_all <- matrix(0,nrow=p,ncol=j+1)
	V_temp <- matrix(0,nrow=n,ncol=p)
    for (k in (1:p)[-imp]){
      #search the maximum, a_opt is the current A regime
      X_temp <- X[,c(imp,k)]
      reg <- lm(Y~X[,c(imp,k)]*A)
      coef <- reg$coef
      coef_all[k,] <- reg$coef[(j+2):(2*(j+1))]
      diff <- cbind(1,X[,c(imp,k)])%*%coef_all[k,]
      g_opt <- as.numeric(diff>=0)
      S_temp[k] <- sum(diff*(g_opt-a_opt[,j]))
	  V_temp[,k] <- Y+diff*(g_opt-A)
    }
    #important variables
    j_opt<-which.max(S_temp)
    imp<-c(imp,j_opt)
    coef_temp<-coef_all[j_opt,]
    regime_coef[j+1,c(1,imp+1)]<-coef_temp
    a_opt[,j+1]<-as.numeric(cbind(1,X[,imp])%*%coef_temp>0)
    V_est[,j+1]<-V_temp[,j_opt]
	S[j+1]<-S_temp[j_opt]
	progress(100*j/step)
  }
  
  sp <- which.min(-log(cumsum(S))+(0:step)*log(n)/n)
  return (list(beta=as.vector(regime_coef[sp, ]), V=V_est[,sp]))
}
