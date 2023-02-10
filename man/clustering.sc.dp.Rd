\name{clustering.sc.dp}
\alias{clustering.sc.dp}
\title{
Optimal Clustering Multidimensional Data with Sequential Constraint by Dynamic Programming
}
\description{
Perform optimal clustering on multidimensional data with sequential constraint (i.e. only subsequent elements of the input may form a cluster).
}
\usage{
clustering.sc.dp(x, k)
}
\arguments{
  \item{x}{a multi-dimensional array containing input data to be clustered}
  \item{k}{the number of clusters}
}
\details{
The '\code{clustering.sc.dp}' algorithm (Szkaliczki, 2016) groups multidimensional data given by \code{x} into \code{k} clusters with sequential constraint by dynamic programming. It generalises the one-dimensional '\code{Ckmeans.1d.dp}' algorithm (Wang and Song, 2011) to multidimensional data. If only subsequent elements of the input data may form a cluster the algorithm guarantees the optimality of clustering -- the sum of squares of within-cluster distances (\var{withinss}) from each element to its corresponding cluster centre (mean) is always the minimum. The sequential constraint is typically required in clustering datastreams or items with time stamps such as video frames, GPS signals of vehicles or movement data of persons etc. The run time of the algorithm is \var{O}( \var{k} \var{d} \var{n}^2) where \var{k}, \var{d} and \var{n} gives the number of clusters, the dimensions of the elements and the number of elements, respectively. 
}
\value{
	An object of class '\code{clustering.sc.dp}' which has a print method and is a list with components:
	\item{cluster}{a vector of cluster indices assigned to each element in \code{x}. Each cluster is indexed by an integer from 1 to \code{k}}
	\item{centers}{a matrix whose rows represent cluster centres}
	\item{withinss}{the within-cluster sum of squares for each cluster}
	\item{size}{a vector of the number of points in each cluster}
}
\references{
Szkaliczki, T. (2016) "clustering.sc.dp: Optimal Clustering with Sequential Constraint by Using Dynamic Programming" <doi: 10.32614/RJ-2016-022>
Wang, H. and Song, M. (2011) "Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming" <doi: 10.32614/RJ-2011-015>
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
k<-2
result<-clustering.sc.dp(x,k)
plot(x, type = 'b', col = result$cluster)
points(result$centers, pch = 24, bg = (1:k))
}
\keyword{optimize}
\keyword{cluster}
\keyword{ts}
