\name{Rpow2}
\alias{Rpow2}
\title{Provides 'repeated power to the 2' as accessible function}
\description{
Provides 'repeated power to the 2' as accessible function
}
\usage{
Rpow2(z)   }
\arguments{
  \item{z}{a scalar or vector of numerical values}
}
\details{
 The function returns $z^2 * log(z)$.
    }                           
\value{
 z^2 * log(z)
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
a<-c(1,4,6)
Rpow2(a)
  
 # returns: [1]  0.00000 22.18071 64.50334



}
