\name{PAL}
\alias{PAL}
\alias{PAL.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Penalized A-learning for optimal dynamic treatment regime
}
\description{
  Selects important variables that are involved in the optimal treatment regime based on 
  penalized A-learning estimating equation. This function can be applied to two-stage 
  studies where treatments are sequentially assigned at two different time points. 
}
\usage{
PAL(formula, data, subset, na.action, IC = c("BIC", "CIC", "VIC"), 
    lambda.list = exp(seq(-3.5, 2, 0.1)), refit = TRUE, control = PAL.control(...), 
	model = TRUE, y = TRUE, a1 = TRUE, x1 = TRUE, a2 = TRUE, x2 = TRUE, ...)
	
PAL.fit(y, x1, x2 = NULL, a1, a2 = NULL, IC = c("BIC", "CIC", "VIC"), 
    lambda.list = exp(seq(-3.5, 2, 0.1)), refit = TRUE, 
    control = PAL.control())	
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
  \item{IC}{Information criterion used in determining the regularization parameter. See 'Details'.
}
  \item{lambda.list}{A list of regularization parameter values. Default is \code{exp(seq(-3.5, 2, 0.1))}.
}
  \item{refit}{After variable selection, should the coefficients be refitted using A-learning
  estimating equation? Default is TRUE.
}
  \item{control}{A list of control argument via \code{\link{PAL.control}}.
}
  \item{model}{A logical value indicating whether \emph{model frame} should be included
  as a component of the return value.
}
  \item{y, a1, x1, a2, x2}{For \code{PAL}: logical values indicating whether the response, 
  the first and second treatments, the baseline and intermediate covariates should be included 
  as a component of the return value.
  
  For \code{PAL.fit}: y is the response vector (the larger the better), a1 and a2 are the first 
  and second treatments patients receive, x1 and x2 are the design matrices consisting of 
  patients' baseline covariates and intermediate covariates.   
}
  \item{\dots}{Argument passed to \code{\link{PAL.control}}.
}
}
\details{Penalized A-learning is developed to select important variables involved in the optimal
individualized treatment regime. An individualized treatment regime is a function that maps patients
covariates to the space of available treatment options. The method can be applied to both single-stage
and two-stage studies.

PAL applied the Dantzig selector on the A-learning estimating equation for variable selection. The 
regularization parameter in the Dantzig selector is chosen according to the information criterion.
Specifically, we provide a Bayesian information criterion (BIC), a concordance information criterion 
(CIC) and a value information criterion (VIC). For illustration of these information criteria, consider
a single-stage study. Assume the data is summarized as \eqn{(Y_i, A_i, X_i), i=1,...,n} where \eqn{Y_i}
is the response of the \eqn{i}-th patient, \eqn{A_i} denotes the treatment that patient receives and 
\eqn{X_i} is the corresponding baseline covariates. Let \eqn{\hat{\pi}_i} and \eqn{\hat{h}_i} denote the
estimated propensity score and baseline mean of the \eqn{i}-th patient. For any linear treatment regime
\eqn{I(x^T \beta>c)}, BIC is defined as
\deqn{BIC=-n\log\left( \sum_{i=1}^n (A_i-\hat{\pi}_i)^2 (Y_i-\hat{h}_i-A_i c-A_i X_i^T \beta)^2 \right)-\|\beta\|_0 \kappa_B,}
where \eqn{\kappa_B=\{\log (n)+\log (p+1) \}/\code{kappa}} and \code{kappa} is the model complexity penalty used in the function \code{PAL.control}. 
VIC is defined as
\deqn{VIC=\sum_{i=1}^n \left(\frac{A_i d_i}{\hat{\pi}_i}+\frac{(1-A_i) (1-d_i)}{1-\hat{\pi}_i} \right)\{Y_i-\hat{h}_i-A_i (X_i^T \beta+c)\}+
\{\hat{h}_i+\max(X_i^T \beta+c,0)\}-\|\beta\|_0 \kappa_V,}
where \eqn{d_i=I(X_i^T \beta>-c)} and \eqn{\kappa_V=n^{1/3} \log^{2/3} (p) \log (\log (n))/\code{kappa}}. 
CIC is defined as
\deqn{CIC=\sum_{i\neq j} \frac{1}{n} \left( \frac{(A_i-\hat{\pi}_i) \{Y_i-\hat{h}_i\} A_j}{\hat{\pi}_i (1-\hat{\pi}_i) \hat{\pi}_j}-
\frac{(A_j-\hat{\pi}_j) \{Y_j-\hat{h}_j\} A_i}{\hat{\pi}_j (1-\hat{\pi}_j) \hat{\pi}_i} \right) I(X_i^T \beta> X_j^T \beta)
-\|\beta\|_0 \kappa_C,} 
where \eqn{\kappa_C=\log (p) \log_{10}(n) \log(\log_{10}(n))/\code{kappa}}. 

Under certain conditions, it can be shown that CIC and VIC is consistent as long as either the estimated
propensity score or the estimated baseline is consistent.

For single-stage study, the formula should specified as y ~ x1 | a1 where y is the reponse vector (y
should be specified in such a way that a larger value of y indicates better clinical outcomes), x1 is
patient's baseline covariates and a1 is the treatment that patient receives. 

For two-stage study, the formula should be specified as y ~ x1 | a1 | x2 | a2 where y is the response
vector, a1 and a2 the vectors of patients' first and second treatments, x1 and x2 are the design matrices
consisting of patients' baseline covariates and intermediate covariates.  

\code{PAL} standardizes the covariates and includes an intercept in the estimated individualized treatment
regime by default. For single-stage study, the estimated treamtent regime is given by \eqn{I(\code{x1}^T \code{beta1.est}>0)}.
For two-stage study, the estimated regime is given by \eqn{\code{a1}=I(x1^T \code{beta1.est}>0)} and \eqn{\code{a2}=I(\code{x}^T \code{beta2.est}>0)}
where \code{x=c(x1, a1, x2)}. 
}
\value{
  \item{beta2.est}{Estimated coefficients in the second decision rule.}
  \item{beta1.est}{Estimated coefficients in the first decision rule.}
  \item{pi2.est}{Estimated propensity score at the second stage.}
  \item{pi1.est}{Estimated propensity score at the first stage.}
  \item{h2.est}{Estimated baseline function at the second stage.}
  \item{h1.est}{Estimated baseline function at the first stage.}
  \item{alpha2.est}{Regression coefficients in the estimated propensity score at the second stage.}
  \item{alpha1.est}{Regression coefficients in the estimated propensity score at the first stage.}
  \item{theta2.est}{Regression coefficients in the estimated baseline function at the second stage.}
  \item{theta1.est}{Regression coefficients in the estimated baseline function at the first stage.}
  \item{model}{The full model frame (if \code{model = TRUE}).}
  \item{y}{Response vector (if \code{y = TRUE}).}
  \item{x1}{Baseline covariates (if \code{x1 = TRUE}).}
  \item{a1}{A vector of first treatment (if \code{a1 = TRUE}).}
  \item{x2}{Intermediate covariates (if \code{x2 = TRUE}).}
  \item{a2}{A vector of second treatment (if \code{a2 = TRUE}).}
}
\references{
Shi, C. and Fan, A. and Song, R. and Lu, W. (2018) High-Dimensional A-Learing for Optimal
Dynamic Treatment Regimes. \emph{Annals of Statistics,} in press. 

Shi, C. and Song, R. and Lu, W. (2018) Concordance and Value Information Criteria for
Optimal Treatment Decision. \emph{Under review}. 
}
\author{
Chengchun Shi and Ailin Fan
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{PAL.control}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\donttest{
## single-stage study
set.seed(12345)
n <- 200
p <- 1000
X <- matrix(rnorm(n*p), nrow=n, ncol=p)
A <- rbinom(n, 1, 0.5)
CX <- (X[,1] + X[,2])
h <- 1 + X[,1] * X[,3]
Y <- h + A*CX + 0.5*rnorm(n)
result <- PAL(Y~X|A)

## two-stage study
set.seed(12345*2)
n <- 200
p <- 1000
X1 <- matrix(rnorm(n*p), nrow=n, ncol=p)
A1 <- rbinom(n, 1, 0.5)
X2 <- X1[,1] + A1 + 0.5*rnorm(n)
A2 <- rbinom(n, 1, 0.5)
Y <- A2*(A1 + X2) + A1*X1[,1] + 0.5*rnorm(n)
result <- PAL(Y~X1|A1|X2|A2)
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
result <- PAL(Y~X|A)

## two-stage study
set.seed(12345*2)
n <- 50
p <- 20
X1 <- matrix(rnorm(n*p), nrow=n, ncol=p)
A1 <- rbinom(n, 1, 0.5)
X2 <- X1[,1] + A1 + 0.5*rnorm(n)
A2 <- rbinom(n, 1, 0.5)
Y <- A2*(A1 + X2) + A1*X1[,1] + 0.5*rnorm(n)
result <- PAL(Y~X1|A1|X2|A2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Penalized A-learning}
