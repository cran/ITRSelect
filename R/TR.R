TR <-
function(object, x1, a1=NULL, x2=NULL, stage=1){
	if (class(object)=="PAL" || class(object)=="SAS"){
		if (stage==1)
			return (as.vector(((cbind(1, x1) %*% object$beta1.est) >0)+0))
		else
			return (as.vector(((cbind(1, x1, a1, x2) %*% object$beta2.est) >0)+0))
	}
	else {
		stop("Object must be fitted by PAL or SAS!")
	}
}
