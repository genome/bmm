\name{bmm.plot.2d}

\alias{bmm.plot.2d}

\title{
   bmm.plot.2d:  Plot data highlighted according to clustering
}

\description{
   Plot the data with each item's respective color/symbol indicating
   its cluster.  Also display "1-sigma" contour intervals for the
   standard error of the mean for the standard error.
}


\usage{
  bmm.plot.2d(X, mu, alpha, nu, beta, E.pi, r, title, xlab, ylab)
}

\arguments{
   \item{X}{
     an N x D matrix with rows being the items to cluster.
     All entries are assumed to be proportions (i.e., between 0 and 1).
     Notice that there are no summation restrictions--i.e., proportions
     do not sum to unity across an item's dimensions.
   }
   \item{mu}{
     a D x N.c matrix holding the values of the 
     shape parameters for the gamma prior distributions over the 
     u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
     Introduced in eqn (15).
   }
   \item{alpha}{
     a D x N.c matrix holding the values of the 
     rate (i.e., inverse scale) parameters for the gamma prior 
     distributions over the u parameters.  i.e., mu[d,n] is the rate
     parameter governing u[d,n]. Introduced in eqn (15).
   }
   \item{nu}{
     a D x N.c matrix holding the values of the 
     shape parameters for the gamma prior distributions over the 
     v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
     Introduced in eqn (16).
   }
   \item{beta}{
     a D x N.c matrix holding the values of the 
     rate (i.e., inverse scale) parameters for the gamma prior 
     distributions over the v parameters.  i.e., beta[d,n] is the rate
     parameter governing v[d,n]. Introduced in eqn (16).
   }
   \item{E.pi}{
     the D-vector holding the values E[pi], i.e., the expected values
     of the mixing coefficients, defined in eqn (53).
   }
   \item{r}{  
     the N x N.c matrix of responsibilities, with r[n, nc] giving the
     probability that item n belongs to component nc
   }
   \item{title}{
     plot title
   }
   \item{xlab}{
     x label
   }
   \item{ylab}{  
     y label
   }
}

\details{
}

\value{
   ggplot object displaying the clustered data and "1-sigma" contour intervals
   for the standard error of the mean and for the standard error
}

\examples{
}


