\name{censNID}
\alias{censNID}

\title{mle NID sample with single censor point}

\description{Uses AS 138}

\usage{censNID(y, n, cy, ctyp = c("left", "right"))}

\arguments{
  \item{y}{complete observations}
  \item{n}{sample size}
  \item{cy}{censor point}
  \item{ctyp}{either "left" or "right", default is "left"}
}


\value{
a list with the following components:
\item{est}{matrix with mean and sd and their se}
\item{covMat}{covariance matrix}
\item{nobs}{number of complete observations}
\item{itercount}{number of iterations}
\item{ifault}{fault indicator. 0-converged. -1:convergence not obtained in default number of allowed iterations. -2:fewer than 2 complete obserations.
-3:error with confined estimator initial values (should not occur in R version). -4:invalid specification for limits with confined observation. }
}

\references{
M. S. Wolynetz (1979). 
Algorithm AS 138: Maximum Likelihood Estimation from Confined and Censored Normal Data. 
Journal of the Royal Statistical Society. Series C (Applied Statistics), 
Vol. 28, No. 2, pp. 185-195
}

\author{
A. I. McLeod
aimcleod@uwo.ca
}

\seealso{
\code{\link{cenNID}}
}

\examples{
#Test data
z <- rep((-2):2, rep(3, 5))
cy <- -1.5
y <- z[z>= cy]
censNID(y, length(z), cy) 

#Test left-censoring
set.seed(32173217)
n <- 100
z <- rnorm(n)
cy <- -1
ind <- z > cy
y <- z[ind]
censNID(y, n, cy)

#Test right-censoring
set.seed(32173217)
n <- 100
z <- rnorm(n)
cy <- 1
ind <- z < cy
y <- z[ind]
censNID(y, n, cy, ctyp="right")


}

\keyword{ htest }
\keyword{ models }



