\name{Rsqrt}
\alias{Rsqrt}
\title{Provides 'repeated square root' as accessible function}
\description{
Provides 'repeated square root' as accessible function
}
\usage{
Rsqrt(z)  }
\arguments{
  \item{z}{a scalar or vector of numerical values}
}
\details{
 The function returns $sqrt(z) * log(z)$.
    }                           
\value{
 sqrt(z) * log(z)
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
a<-c(1,4,6)
Rsqrt(a)
 # returns: [1] 0.000000 2.772589 4.388896
}
