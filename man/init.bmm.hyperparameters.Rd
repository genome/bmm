\name{init.bmm.hyperparameters}

\alias{init.bmm.hyperparameters}

\title{
  init.bmm.hyperparameters:  Provides default initialization of hyperparameters for Beta mixture model.
}

\description{
  Initialize the bmm hyperparameters (to be passed to bmm or 
  bmm.fixed.num.components)

  NB:  This provides what should be a generally reasonable initialization of
       hyperparameters.  However, better results may be obtained by
       tuning these in an application-specific manner.
}

\usage{
  init.bmm.hyperparameters(X, N.c)
}

\arguments{
  \item{X}{
       an N x D matrix with rows being the items to cluster.
       All entries are assumed to be proportions (i.e., between 0 and 1).
       Notice that there are no summation restrictions--i.e., proportions
       do not sum to unity across an item's dimensions.
  }
  \item{N.c}{ 
       N.c:  the number of components/clusters to attempt
  }
}

\details{
}

\value{
  Returns a list with the following components:
  \item{mu0}{  
      a D x N.c matrix holding the hyperparameter values of the 
      shape parameters for the gamma prior distributions over the 
      u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
      Introduced in eqn (15).
  }
  \item{alpha0}{ 
      a D x N.c matrix holding the hyperparameter values of the 
      rate (i.e., inverse scale) parameters for the gamma prior 
      distributions over the u parameters.  i.e., mu[d,n] is the rate
      parameter governing u[d,n]. Introduced in eqn (15).
  }
  \item{nu0}{ 
      a D x N.c matrix holding the hyperparameter values of the 
      shape parameters for the gamma prior distributions over the 
      v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
      Introduced in eqn (16).
  }
  \item{beta0}{ 
      a D x N.c matrix holding the hyperparameter values of the 
      rate (i.e., inverse scale) parameters for the gamma prior 
      distributions over the v parameters.  i.e., beta[d,n] is the rate
      parameter governing v[d,n]. Introduced in eqn (16).
  }
  \item{c0}{ 
      a vector with D components holding the hyperparameter values of the
      parameters of the Dirichlet distribution over the mixing 
      coefficients pi.  Introduced in eqn (19).
  }
}

\examples{
}


