\name{cenNID}
\alias{cenNID}

\title{mle estimation for censored normal}

\description{
MLE estimation using AS 138 for censored and confined normal.
}

\usage{cenNID(y, L=rep(-Inf, length(y)), U=rep(Inf,length(y))) }

\arguments{
  \item{y}{observations}
  \item{L}{lower bound for each observation that is equal to -Inf when none, default -Inf}
  \item{U}{upper bound for each observation that is equal to Inf when none, default Inf }
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
\code{\link{censNID}}
}

\examples{
set.seed(32173217)
n<-100
z <- rnorm(n)
U <- rep(Inf, n)
L <- rep(-Inf, n)
cy <- -1
ind <- z > cy
m <- sum(as.integer(ind))
y <- c(z[ind], rep(cy, n-m))
L[(m+1):n] <- cy
cenNID(y, L, U) 
}

\keyword{ htest }
\keyword{ models }


