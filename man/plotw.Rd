\name{plotw}
\alias{plotw}
\title{Plot weights of weighted estimation in Cox regression}
\description{
             This function plots the weights used in a weighted Cox regression analysis against time.
}
\usage{
plotw(object, rank=F, log=F, ...)}
\arguments{
  \item{object}{a coxphw object}
  \item{rank}{if set to TRUE, plots the weights against ranked time (default=F)}
  \item{log}{if set to TRUE, shows logarithm of weights (default=F)}
  \item{...}{additional parameters}
}
\details{
 The function plots the raw, inverse-censoring-Kaplan-Meier, and normalized total weights against (ranked) time.
    }
\value{
}
\author{Georg Heinze}
\seealso{coxphw}
\keyword{survival}
