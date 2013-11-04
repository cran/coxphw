\name{PT}
\alias{PT}
\title{Pretransformation function}
\description{
Provides automatic pretransformation of variables (to well-scaled an nonzero values)
}
\usage{
PT(z)      }
\arguments{
  \item{z}{a vector of numerical values}
}
\details{
 The function transforms a variable by shifting to positive values, and dividing by scaling factor (a power of 10) such that the SD is approximately equal to 1.
    }
\value{
 (a+shift)/scale
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
 a<-c(-6,-1,4,6)
 PT(a)
 # returns: [1] 0.2 0.7 1.2 1.4
 
}
