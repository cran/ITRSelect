\name{PAL.control}
\alias{PAL.control}
\title{Control parameters for penalized A-learning}
\description{
  Parameters that control fitting of penalized A-learning.
}
\usage{
PAL.control(pi1.est = NULL, pi2.est = NULL, h1.est = NULL, h2.est = NULL, kappa = NULL, 
  penalty = 'SCAD')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pi1.est}{Estimated propentisy score at the first stage. By default, a penalized 
  logistic regression model is fitted to estimate the propensity score.
}
  \item{pi2.est}{Estimated propentisy score at the second stage. By default, a penalized 
  logistic regression model is fitted to estimate the propensity score.
}
  \item{h1.est}{Estimated baseline function at the first stage. By default, a penalized
  linear regression model is fitted to estimate the baseline function.
}
  \item{h2.est}{Estimated baseline function at the second stage. By default, a penalized 
  linear regression model is fitted to estimate the baseline function.
}
  \item{kappa}{The model complexity penalty used in the information criteria. By default, 
  \eqn{kappa=1} if BIC or CIC is used and \eqn{kappa=4} if VIC is used. 
}
  \item{penalty}{The penalty to be applied to the propensity score and baseline model.
  Either "MCP", "SCAD" (the default), or "lasso". 
}
}
\value{
  A list with the arguments specified. 
}

\seealso{\code{\link{PAL}}, \code{\link{PAL.fit}}
}
\examples{
\donttest{
set.seed(12345)
n <- 200
p <- 200
X <- matrix(rnorm(n*p), nrow=n, ncol=p)
A <- rbinom(n, 1, 0.5)
CX <- (X[,1] + X[,2])
h <- 1 + X[,1] * X[,3]
Y <- h + A*CX + 0.5*rnorm(n)
result <- SAS(Y~X|A, pi1.est=0.5)
}
}

