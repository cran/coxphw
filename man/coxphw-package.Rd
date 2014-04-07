\name{coxphw-package}

\alias{coxphw-package}

\docType{package}

\title{ Weighted Estimation in Cox Regression }

\description{This package implements weighted estimation in Cox regression as proposed by 
             Schemper, Wakounig and Heinze (Statistics in Medicine, 2009). Weighted Cox regression 
             provides unbiased average hazard ratio estimates also in case of non-proportional hazards. 
             Additionally estimation of nonlinear effects using fractional polynomials similar to the 
             MFP algorithm (Royston, Sauerbrei, 2008) is provided. This feature can also be used to 
             estimate the interaction of a covariate with a nonlinear function of time. }

\details{
\tabular{ll}{
Package: \tab \code{coxphw}\cr
Type: \tab Package\cr
Version: \tab 3.0.0\cr
Date: \tab 2014-04-07\cr
License: \tab GPL-2\cr
}

Main functions included in the \code{coxphw} package are
\tabular{ll}{
\code{\link{coxphw}}    \tab weigthed estimation of Cox regression: either (recommended) estimation of \cr 
                        \tab average hazard ratios (Schemper et al., 2009), estimation of average regression \cr
                        \tab effects (Xu and O'Quigley, 2000), or proportional hazards regression. \cr
                        \tab Estimation of nonlinear effects using fractional polynomials is possible. \cr 
\code{\link{concord}}   \tab obtain generalized concordance probabilities with confidence intervalls. \cr 
\code{\link{plotw}}     \tab plots the weights used in a weighted Cox regression analysis against time. \cr 
\code{\link{fp.power}}  \tab provides fractional polynomials as accessible function. \cr 
\code{\link{plotshape}} \tab either visualizes a nonlinear or a time-dependent effect (possibly estimated \cr
                        \tab by fractional polynomials) of a \code{coxphw} fit or computes the effect estimates \cr
                        \tab at pre-specified values. It plots the relative or log relative hazard versus \cr
                        \tab values of a continuous covariable. \cr 
\code{\link{wald}}      \tab obtain Wald chi-squared test statistics and p-values for one or more regression \cr
                        \tab coefficients given their variance-covariance matrix. \cr                           
}

Data sets included in the \code{coxphw} package are
\tabular{ll}{
\code{biofeedback} \tab biofeedback treatment data \cr 
\code{gastric}     \tab gastric cancer data \cr 
}
}

\note{ The \code{SAS} macro \code{WCM} with similar functionality (but without fractional polynomials) 
       can be obtained at \cr
       \url{http://cemsiis.meduniwien.ac.at/en/kb/science-research/software/statistical-software/wcmcoxphw/}. 
       
       Up to Version 2.13 \pkg{coxphw} used a slightly different syntax (arguments: \code{AHR}, \code{AHR.norobust}, 
       \code{ARE}, \code{PH}, \code{normalize}, \code{censcorr}, \code{prentice}, \code{breslow}, 
       \code{taroneware}). From Version 3.0 on the old syntax is disabled. }

\author{R by Meinhard Ploner, Georg Heinze, Daniela Dunkler; Fortran by Georg Heinze \cr
        {Maintainer: \email{georg.heinze@meduniwien.ac.at}} }

\references{Dunkler D, Schemper M, Heinze G. (2010) Gene Selection in Microarray Survival Studies Under 
Possibly Non-Proportional Hazards. \emph{Bioinformatics} \bold{26}:784-90.

Lin D and Wei L (1989). The Robust Inference for the Cox Proportional Hazards Model. 
\emph{Journal of the American Statistical Association} \bold{84}, 1074-1078.

Lin D (1991). Goodness-of-Fit Analysis for the Cox Regression Model Based on a Class of Parameter 
Estimators. \emph{Journal of the American Statistical Association} \bold{86}, 725-728.

Royston P and Altman D (1994). Regression Using Fractional Polynomials of Continuous Covariates: 
Parsimonious Parametric Modelling. \emph{Applied Statistics} \bold{43}, 429-467.

Royston P and Sauerbrei W (2008). \emph{Multivariable Model-Building. A Pragmatic Approach to Regression 
Analysis Based on Fractional Polynomials for Modelling Continuous Variables.} Wiley, Chichester, UK.

Sasieni P (1993). Maximum Weighted Partial Likelihood Estimators for the Cox Model. \emph{Journal of 
the American Statistical Association} \bold{88}, 144-152.

Schemper M (1992).  Cox Analysis of Survival Data with Non-Proportional Hazard Functions. \emph{The 
Statistician} \bold{41}, 455-465.

Schemper M, Wakounig S and Heinze G (2009). The Estimation of Average Hazard Ratios by Weighted Cox 
Regression. \emph{Statistics in Medicine} \bold{28}, 2473-2489.

Xu R and O'Quigley J (2000). Estimating Average Regression Effect Under Non-Proportional Hazards. 
\emph{Biostatistics} \bold{1}, 423-439. }

\seealso{ \code{\link{coxphw}}, \code{\link{concord}}, \code{\link{plotshape}}, \code{\link{plotw}}, \code{\link{fp.power}}, \code{\link{wald}} }

\examples{
## for examples see coxphw
}

\keyword{ package }
\keyword{ regression }
\keyword{ survival }
