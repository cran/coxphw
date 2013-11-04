\name{powM0.5}
\alias{powM0.5}
\title{Provides 'power to the minus 0.5' as accessible function}
\description{
Provides 'power to the minus 0.5' as accessible function
}
\usage{
powM0.5(z)    }
\arguments{
  \item{z}{a scalar or vector of numerical values}
}
\details{
 The function returns $z^{-0.5}$.
    }                           
\value{
 z^{-0.5}
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
 a<-c(1,4,6)
 powM0.5(a)

 # returns: [1] 1.0000000 0.5000000 0.4082483

}
