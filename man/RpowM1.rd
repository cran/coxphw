\name{RpowM1}
\alias{RpowM1}
\title{Provides 'repeated power to the minus 1' as accessible function}
\description{
Provides 'repeated power to the minus 1' as accessible function
}
\usage{
RpowM1(z)  }
\arguments{
  \item{z}{a scalar or vector of numerical values}
}
\details{
 The function returns $z^(-1) * log(z)$.
    }                           
\value{
 z^(-1) * log(z)
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
a<-c(1,4,6)
RpowM1(a)
  
 # returns: [1] 0.0000000 0.3465736 0.2986266



}
