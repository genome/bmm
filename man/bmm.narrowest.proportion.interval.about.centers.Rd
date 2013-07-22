\name{bmm.narrowest.proportion.interval.about.centers}

\alias{bmm.narrowest.proportion.interval.about.centers}

\title{
  bmm.narrowest.proportion.interval.about.centers:  Compute intervals on proportions
}

\description{
  bmm.narrowest.proportion.interval.about.centers:  Return the narrowest 
    interval of proportions from the posterior distribution that
    includes their centers a fixed proportion of the posterior distribution
    (e.g., .68 or .95).
}

\usage{
  bmm.narrowest.proportion.interval.about.centers(mu, alpha, nu, beta, proportion) 

}

\arguments{
   \item{mu}{
        a D x N.c matrix holding the shape parameters for the gamma prior 
        distributions over the u parameters.  i.e., mu[d,n] is the 
        shape parameter governing u[d,n].
        Introduced in eqn (15).
   }
   \item{alpha}{
        a D x N.c matrix holding the rate (i.e., inverse scale) 
        parameters for the gamma prior distributions over the u 
        parameters.  i.e., mu[d,n] is the rate parameter governing u[d,n].
        Introduced in eqn (15).
   }
   \item{nu}{
        a D x N.c matrix holding the hyperparameter values of the 
        shape parameters for the gamma prior distributions over the 
        v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
        Introduced in eqn (16).
   }
   \item{beta}{
        a D x N.c matrix holding the hyperparameter values of the 
        rate (i.e., inverse scale) parameters for the gamma prior 
        distributions over the v parameters.  i.e., beta[d,n] is the rate
        parameter governing v[d,n]. Introduced in eqn (16).
   }
   \item{proportion}{
        the percentage of the distribution of each proportion to include
        in the interval (e.g., .68 or .95)
   }
}

\details{
}

\value{
  a list with values lb and ub, where each are D x N.c
  matrices holding the respective lower or upper bounds of the intervals,
  and centers, which is a D x N.c matrix holding the empirical means sampled
  from the posterior distribution
}

\examples{
}




