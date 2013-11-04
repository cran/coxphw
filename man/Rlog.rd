\name{Rlog}
\alias{Rlog}
\title{Provides 'repeated logarithm' as accessible function}
\description{
Provides 'repeated logarithm' as accessible function
}
\usage{
Rlog(z)   }
\arguments{
  \item{z}{a scalar or vector of numerical values}
}
\details{
 The function returns $log(z) * log(z)$.
    }                           
\value{
 log(z) * log(z)
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
a<-c(1,4,6)
Rlog(a)
  
 # returns: [1] 0.000000 1.921812 3.210402
}
