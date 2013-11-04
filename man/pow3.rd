\name{pow3}
\alias{pow3}
\title{Provides 'power to the 3' as accessible function}
\description{
Provides 'power to the 3' as accessible function
}
\usage{
pow3(z)  }
\arguments{
  \item{z}{a scalar or vector of numerical values}
}
\details{
 The function returns $z^3$.
    }                           
\value{
 z^3
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
 a<-c(1,4,6)
 pow3(a)

 # returns: [1]   1  64 216

}
