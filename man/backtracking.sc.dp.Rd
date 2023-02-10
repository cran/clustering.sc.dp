\name{backtracking.sc.dp}
\alias{backtracking.sc.dp}
\title{
Backtracking Clustering for a Specific Cluster Number
}
\description{
Creates clustering for \code{k} number of clusters by using the backtrack data produced by \code{findwithinss.sc.dp()}. 
}
\usage{
backtracking.sc.dp(x, k, backtrack)
}
\arguments{
  \item{x}{ a multi-dimensional array containing input data to be clustered}
  \item{k}{ the number of clusters}
  \item{backtrack}{the backtrack data}
}
\details{
    If the number of clusters is unknown \code{findwithinss.sc.dp()} followed by \code{backtracking.sc.dp()} can be used for performing clustering. If only subsequent elements of the input data may form a cluster method \code{findwithinss.sc.dp()} calculates the exact minimum of the sum of squares of within-cluster distances (\var{withinss}) from each element to its corresponding cluster centre (mean) for different cluster numbers. The user may analyse the \var{withinss} in order to select the proper number of clusters. In this case, it is enough to run method \code{backtracking.sc.dp()} only once. Another option is to run \code{findwithinss.sc.dp()} once, repeat the \code{backtracking.sc.dp()} step for a range of potential cluster numbers and then the user may evaluate the optimal solutions created for different number of clusters. This requires much less time than repeating the whole clustering algorithm for the different cluster numbers. 
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

\seealso{
\code{\link{findwithinss.sc.dp}}, \code{\link{clustering.sc.dp}}
}
\examples{
# Example: clustering data generated from a random walk with small withinss
x<-matrix(, nrow = 100, ncol = 2)
x[1,]<-c(0,0)
for(i in 2:100) {
  x[i,1]<-x[i-1,1] + rnorm(1,0,0.1)
  x[i,2]<-x[i-1,2] + rnorm(1,0,0.1)
}
k<-10
r<-findwithinss.sc.dp(x,k)

# select the first cluster number where withinss drops below a threshold
thres <- 5.0
k_th <- 1;
while(r$twithinss[k_th] > thres & k_th < k) {
    k_th <- k_th + 1
}

# backtrack
result<-backtracking.sc.dp(x,k_th, r$backtrack)
plot(x, type = 'b', col = result$cluster)
points(result$centers, pch = 24, bg = (1:k_th))
}
\keyword{optimize}
\keyword{cluster}
\keyword{ts}
