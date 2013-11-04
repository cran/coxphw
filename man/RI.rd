\name{RI}
\alias{RI}
\title{Provides 'repeated identity' as accessible function}
\description{
Provides 'repeated identity' as accessible function
}
\usage{
RI(z)   }
\arguments{
  \item{z}{a scalar or vector of numerical values}
}
\details{
 The function returns $z * log(z)$.
    }                           
\value{
 z * log(z)
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
a<-c(1,4,6)
 RI(a)
  
 # returns: [1]  0.000000  5.545177 10.750557


}
