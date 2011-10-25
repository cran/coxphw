\name{RpowM2}
\alias{RpowM2}
\title{Provides 'repeated power to the minus 2' as accessible function}
\description{
Provides 'repeated power to the minus 2' as accessible function
}
\usage{
RpowM2(z)      }
\arguments{
  \item{z}{a scalar or vector of numerical values}
}
\details{
 The function returns $z^(-2) * log(z)$.
    }                           
\value{
 z^(-2) * log(z)
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
a<-c(1,4,6)
RpowM2(a)
  
 # returns: [1] 0.0000000 0.0866434 0.0497711


}
