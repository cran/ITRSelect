\name{TR}
\alias{TR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Individualized treatment regime based on PAL or SAS.
}
\description{
Recommend individualized treatment regime for future patients, based on the 
penalized A-learning method, or sequential advantage selection method.
}
\usage{
TR(object, x1, a1 = NULL, x2 = NULL, stage = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Fitted object of class "PAL" or "SAS".
}
  \item{x1}{A matrix consisting of future patients baseline covariates.
}
  \item{a1}{A vector consisting of future patients first treatments. Not
  needed if \code{stage = 2}.
}
  \item{x2}{A matrix consisting of future patients intermediate covariates.
  Not needed if \code{stage = 2}. 
}
  \item{stage}{Outputs the first-stage decision rule for future patients if \code{stage = 1}.
  Otherwise, outputs the second-stage decision rule for future patients.
  }
}
\value{
A vector of individualized treatments tailored for future patients.
}
\author{
Chengchun Shi
}
\seealso{\code{\link{PAL}}, \code{\link{SAS}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\donttest{
## load simulated STARD data
data(SSTARD.twostage)
## estimate individualized treatment regime using SAS
result <- SAS(Y~X1|A1|X2|A2, data=SSTARD.twostage)
## make the recommendation 
TR(result, x1=SSTARD.twostage$X1, a1=SSTARD.twostage$A1, x2=SSTARD.twostage$X2, stage=2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Individualized treatment regime}
