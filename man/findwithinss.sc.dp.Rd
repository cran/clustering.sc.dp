\name{findwithinss.sc.dp}
\alias{findwithinss.sc.dp}
\title{
Finding Optimal Withinss in Clustering Multidimensional Data with Sequential Constraint by Dynamic Programming
}
\description{
Performs the main step of clustering multidimensional data with sequential constraint by a dynamic programming approach guaranteeing optimality. 
It returns the minimum \var{withinss} for each number of clusters less than or equal to \code{k} and backtracking data that can be used to find quickly the optimal clustering for a specific cluster number. This function was created in order to support the case when the number of clusters is unknown in advance. 
}
\usage{
findwithinss.sc.dp(x, k)
}
\arguments{
  \item{x}{	a multi-dimensional array containing input data to be clustered}
  \item{k}{	the maximal number of clusters, the output will be generated for cluster numbers between 1 and \code{k} } 
}
\details{
	Method \code{clustering.sc.dp()} is split into two methods (\code{findwithinss.sc.dp()} and \code{backtracking.sc.dp()}) in order to support the case when the number of clusters is not known in advance. Method \code{findwithinss.sc.dp()} returns the minimal sum of squares of within-cluster distances (\var{withinss}) for each number of clusters less than or equal to \code{k} and the backtrack data which can be used to quickly determine the optimal clustering for a specific cluster number. The returned \var{withinss} are guaranteed to be optimal among the solutions where only subsequent items form a cluster. 

	The outputs of the method can be used to select the proper number of clusters. The user may analyse the \var{withinss} in order to select the proper number of clusters. Another option is to run \code{findwithinss.sc.dp()} once, repeat the \code{backtracking.sc.dp()} step for a range of potential cluster numbers and then the user may evaluate the optimal solutions created for different number of clusters. This requires much less time than repeating the whole clustering algorithm. 
	}
\value{
	A list with components:
	\item{twithinss}{ a vector of total within-cluster sums of the optimal clusterings for each number of clusters less than or equal to \code{k}.}
	\item{backtrack}{ backtrack data used by \code{backtracking.sc.dp()}.}
}
\author{
	Tibor Szkaliczki \email{szkaliczki.tibor@sztaki.hu}
}

\seealso{
\code{\link{clustering.sc.dp}}, \code{\link{backtracking.sc.dp}}
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
