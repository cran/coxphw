\name{Rpow3}
\alias{Rpow3}
\title{Provides 'repeated power to the 3' as accessible function}
\description{
Provides 'repeated power to the 3' as accessible function
}
\usage{
Rpow3(z)  }
\arguments{
  \item{z}{a scalar or vector of numerical values}
}
\details{
 The function returns $z^3 * log(z)$.
    }                           
\value{
 z^3 * log(z)
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
a<-c(1,4,6)
Rpow3(a)
  
 # returns: [1]   0.00000  88.72284 387.02005



}
