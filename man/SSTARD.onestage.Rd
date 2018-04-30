\name{SSTARD.onestage}
\alias{SSTARD.onestage}
\docType{data}
\title{Simulated single-stage dataset from the STAR*D study
}
\description{
This data is constructed based on the Sequenced Treatment Alternatives to Relieve Depression
(STAR*D) study. Due to a data confidentiality agreement, we are not able to provide the 
original dataset. Here, we generate simulated data that are similar to the original dataset. 
}
\usage{data("SSTARD.onestage")}
\format{
  A list with 319 observations.
  \describe{
	\item{\code{Y}}{Patients' responses. The larger the better.}
	\item{\code{A}}{Treatments patients receive.}
	\item{\code{X}}{A 319*305 matrix consisting of patients' baseline covariates.}
  }
}
%%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
\references{
Fava, M., Rush, A. J., Trivedi, M. H., Nierenberg, A. A., Thase, M. E., Sack-
eim, H. A., Quitkin, F. M., Wisniewski, S., Lavori, P. W., Rosenbaum, J. F.
et al. (2003). Background and rationale for the Sequenced Treatment Alternatives to
Relieve Depression (STAR* D) study. \emph{Psychiatric Clinics of North America,} 
\bold{ :26} 457-494.

Rush, A. J., Fava, M., Wisniewski, S. R., Lavori, P. W., Trivedi, M. H., Sack-
eim, H. A., Thase, M. E., Nierenberg, A. A., Quitkin, F. M., Kashner, T. M. et al. 
(2004). Sequenced treatment alternatives to relieve depression (STAR* D): rationale 
and design. \emph{Controlled clinical trials,} \bold{ :25} 119-142.
}
\examples{
data(SSTARD.onestage)
}
\keyword{datasets}
