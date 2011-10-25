\name{pow2}
\alias{pow2}
\title{Provides 'power to the 2' as accessible function}
\description{
Provides 'power to the 2' as accessible function
}
\usage{
pow2(z)  }
\arguments{
  \item{z}{a scalar or vector of numerical values}
}
\details{
 The function returns $z^2$.
    }                           
\value{
 z^2
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
 a<-c(1,4,6)
 pow2(a)

 # returns: [1]  1 16 36

}
