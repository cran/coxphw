\name{fp.power}
\alias{fp.power}
\title{Provides fractional polynomials as accessible function}
\description{
Provides fractional polynomials as accessible function
}
\usage{
fp.power(z,a,b=NULL)  }
\arguments{
  \item{z}{a scalar or vector of positive numerical values}
  \item{a}{first power}
  \item{b}{optional second power}
}
\details{
 The function returns fp(a) of z (and optionally fp(b) of z).
    }                           
\value{
 a matrix with one or two columns (if a second power b was specified), and number of rows equal to the length of z.
 The columns are sorted by descending power.
}
\author{Georg Heinze}
\seealso{coxph}
\keyword{survival}
\examples{
 z<-c(1,4,6)
 fp.power(z,1)
 fp.power(z,0.5)
 fp.power(z,0.5,0.5)
 fp.power(z,0,2)

}
