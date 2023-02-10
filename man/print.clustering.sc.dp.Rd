\name{print.clustering.sc.dp}
\alias{print.clustering.sc.dp}
\title{
Print the result returned by calling clustering.sc.dp
}
\description{
Print the result returned by calling clustering.sc.dp
}
\usage{
\method{print}{clustering.sc.dp}(x, ...)
}
\arguments{
  \item{x}{ object returned by calling clustering.sc.dp()}
  \item{\dots}{	ignored arguments }
}

\value{
	An object of class '\code{clustering.sc.dp}' which has a print method and is a list with components:
	\item{cluster}{ A vector of integers (\code{1:k}) indicating the cluster to which each point is allocated.}
	\item{centers}{ A matrix whose rows represent cluster centres.}
	\item{withinss}{ The within-cluster sum of squares for each cluster.}
	\item{size}{ The number of points in each cluster.}
}
\author{
	Tibor Szkaliczki \email{szkaliczki.tibor@sztaki.hu}
}
\examples{
# Example: clustering data generated from a random walk
x<-matrix(, nrow = 100, ncol = 2)
x[1,]<-c(0,0)
for(i in 2:100) {
  x[i,1]<-x[i-1,1] + rnorm(1,0,0.1)
  x[i,2]<-x[i-1,2] + rnorm(1,0,0.1)
}
result<-clustering.sc.dp(x,2)
print(result)
}

\keyword{optimize}
\keyword{cluster}
\keyword{ts}
