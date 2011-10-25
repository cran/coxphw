\name{powM1}
\alias{powM1}
\title{Provides 'power to the minus 1' as accessible function}
\description{
Provides 'power to the minus 1' as accessible function
}
\usage{
powM1(z)   }
\arguments{
  \item{z}{a scalar or vector of numerical values}
}
\details{
 The function returns $z^{-1}$.
    }                           
\value{
 z^{-1}
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
 a<-c(1,4,6)
 powM1(a)

 # returns: [1] 1.0000000 0.2500000 0.1666667
 
}
