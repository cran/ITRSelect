\name{SAS}
\alias{SAS}
\alias{SAS.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Sequential advantage selection for optimal dynamic treatment regime
}
\description{
  Select variables that are qualitatively interacted with the treatment based on
  a modified S-score method and a BIC-type criterion. This function can be applied to two-stage 
  studies where treatments are sequentially assigned at two different time points. 
}
\usage{
SAS(formula, data, subset, na.action, step,
	model = TRUE, y = TRUE, a1 = TRUE, x1 = TRUE, a2 = TRUE, x2 = TRUE, ...)
	
SAS.fit(y, x1, x2 = NULL, a1, a2 = NULL, step)	
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A symbolic description of the model to be fitted(of type y ~ x1 | a1 or
  y ~ x1 | a1 | x2 | a2. Details are given 'Details').
}
  \item{data}{An optional list or environment containing variables in \code{formula}.
}
  \item{subset, na.action}{Arguments controlling formula processing via \code{\link{model.frame}}.
}
  \item{step}{SAS uses a forward selection procedure. The maximum size of the model is specified by \code{step}.
  By default, it is equal to \eqn{n/\log(n)} where \eqn{n} is the sample size.
}
  \item{model}{A logical value indicating whether \emph{model frame} should be included
  as a component of the return value.
}
  \item{y, a1, x1, a2, x2}{For \code{SAS}: logical values indicating whether the response, 
  the first and second treatments, the baseline and intermediate covariates should be included 
  as a component of the return value.
  
  For \code{SAS.fit}: y is the response vector (the larger the better), a1 and a2 are the first 
  and second treatments patients receive, x1 and x2 are the design matrices consisting of 
  patients' baseline covariates and intermediate covariates.   
}
  \item{\dots}{Currently not used}.
}
\details{For single-stage study, the formula should specified as y ~ x1 | a1 where y is the reponse vector (y
should be specified in such a way that a larger value of y indicates better clinical outcomes), x1 is
patient's baseline covariates and a1 is the treatment that patient receives. 

For two-stage study, the formula should be specified as y ~ x1 | a1 | x2 | a2 where y is the response
vector, a1 and a2 the vectors of patients' first and second treatments, x1 and x2 are the design matrices
consisting of patients' baseline covariates and intermediate covariates.  

The function returns linear dynamic treatment regimes. For single-stage study, the estimated treamtent regime 
for future patients is given by \eqn{I(\code{x1}^T \code{beta1.est}>0)}. For two-stage study, the estimated regime 
is given by \eqn{\code{a1}=I(x1^T \code{beta1.est}>0)} and \eqn{\code{a2}=I(\code{x}^T \code{beta2.est}>0)}
where \code{x=c(x1, a1, x2)}. 
}
\value{
  \item{beta2.est}{Estimated coefficients in the second decision rule.}
  \item{beta1.est}{Estimated coefficients in the first decision rule.}
  \item{model}{The full model frame (if \code{model = TRUE}).}
  \item{y}{Response vector (if \code{y = TRUE}).}
  \item{x1}{Baseline covariates (if \code{x1 = TRUE}).}
  \item{a1}{A vector of first treatment (if \code{a1 = TRUE}).}
  \item{x2}{Intermediate covariates (if \code{x2 = TRUE}).}
  \item{a2}{A vector of second treatment (if \code{a2 = TRUE}).}
}
\references{
Fan, A. and Lu, W. and Song, R. (2016) Sequential Advantage Selection for Optimal Treatment 
Regime. \emph{Annals of Applied Statistics,} \bold{ 10:} 32-53. 
}
\author{
Ailin Fan and Chengchun Shi
}

\examples{
\donttest{
## single-stage study
set.seed(12345)
n <- 200
p <- 200
X <- matrix(rnorm(n*p), nrow=n, ncol=p)
A <- rbinom(n, 1, 0.5)
CX <- (X[,1] + X[,2])
h <- 1 + X[,1] * X[,3]
Y <- h + A*CX + 0.5*rnorm(n)
result <- SAS(Y~X|A)

## two-stage study
set.seed(12345*2)
n <- 200
p <- 200
X1 <- matrix(rnorm(n*p), nrow=n, ncol=p)
A1 <- rbinom(n, 1, 0.5)
X2 <- X1[,1] + A1 + 0.5*rnorm(n)
A2 <- rbinom(n, 1, 0.5)
Y <- A2*(A1 + X2) + A1*X1[,1] + 0.5*rnorm(n)
result <- SAS(Y~X1|A1|X2|A2)
}

## single-stage study
set.seed(12345)
n <- 50
p <- 20
X <- matrix(rnorm(n*p), nrow=n, ncol=p)
A <- rbinom(n, 1, 0.5)
CX <- (X[,1] + X[,2])
h <- 1 + X[,1] * X[,3]
Y <- h + A*CX + 0.5*rnorm(n)
result <- SAS(Y~X|A)

## two-stage study
set.seed(12345*2)
n <- 50
p <- 20
X1 <- matrix(rnorm(n*p), nrow=n, ncol=p)
A1 <- rbinom(n, 1, 0.5)
X2 <- X1[,1] + A1 + 0.5*rnorm(n)
A2 <- rbinom(n, 1, 0.5)
Y <- A2*(A1 + X2) + A1*X1[,1] + 0.5*rnorm(n)
result <- SAS(Y~X1|A1|X2|A2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Sequential advantage selection}
