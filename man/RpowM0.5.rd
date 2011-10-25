\name{RpowM0.5}
\alias{RpowM0.5}
\title{Provides 'repeated power to the minus 0.5' as accessible function}
\description{
Provides 'repeated power to the minus 0.5' as accessible function
}
\usage{
RpowM0.5(z)   }
\arguments{
  \item{z}{a scalar or vector of numerical values}
}
\details{
 The function returns $z^(-0.5) * log(z)$.
    }                           
\value{
 z^(-0.5) * log(z)
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
a<-c(1,4,6)
RpowM0.5(a)
  
 # returns: [1] 0.0000000 0.6931472 0.7314827
}
