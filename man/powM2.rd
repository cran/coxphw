\name{powM2}
\alias{powM2}
\title{Provides 'power to the minus 2' as accessible function}
\description{
Provides 'power to the minus 2' as accessible function
}
\usage{
powM2(z)      }
\arguments{
  \item{z}{a scalar or vector of numerical values}
}
\details{
 The function returns $z^{-2}$.
    }                           
\value{
 z^{-2}
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
 a<-c(1,4,6)
 powM2(a)
 a
 # returns: [1] 1.00000000 0.06250000 0.02777778
 
}
