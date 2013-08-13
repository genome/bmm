\name{init.bmm.parameters}
\alias{init.bmm.parameters}
\title{
  init.bmm.parameters:  Initialize bmm parameters
}
\description{

   Initialize the bmm parameters (to be passed to bmm or bmm.fixed.num.components

   Initialize parameters such that expected proportions have the values
   determined by an initial k-means clustering.
   
   NB:  This provides what should be a generally reasonable initialization of
        hyperparameters.  However, better results may be obtained by
        tuning these in an application-specific manner.
}

\usage{
   init.bmm.parameters(X, N.c, mu0, alpha0, nu0, beta0, c0, verbose=0)
}
\arguments{
  \item{X}{
       an N x D matrix with rows being the items to cluster.
       All entries are assumed to be proportions (i.e., between 0 and 1).
       Notice that there are no summation restrictions--i.e., proportions
       do not sum to unity across an item's dimensions.
  }
  \item{N.c}{
       the number of components/clusters to attempt
  }
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
  \item{verbose}{
       output progress in terms of mixing coefficient (expected) values if 1.
  }
    
}

\value{
  A list with the following entries:

  \item{mu}{
       a D x N.c matrix holding the _initial_ values of the 
       shape parameters for the gamma prior distributions over the 
       u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
       NB:  this is the initial value mu, which is updated upon iteration.
       It is not (necessarily) the same as the hyperparameter mu0, which
       is unchanged by iteration.
       Introduced in eqn (15).
  }
  \item{alpha}{
       a D x N.c matrix holding the _initial_ values of the 
       rate (i.e., inverse scale) parameters for the gamma prior 
       distributions over the u parameters.  i.e., mu[d,n] is the rate
       parameter governing u[d,n]. Introduced in eqn (15).
       NB:  this is the initial value alpha, which is updated upon iteration.
       It is not (necessarily) the same as the hyperparameter alpha0, which
       is unchanged by iteration.
  }
  \item{nu}{
       a D x N.c matrix holding the _initial_ values of the 
       shape parameters for the gamma prior distributions over the 
       v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
       Introduced in eqn (16).
       NB:  this is the initial value nu, which is updated upon iteration.
       It is not (necessarily) the same as the hyperparameter nu0, which
       is unchanged by iteration.
  }
  \item{beta}{
       a D x N.c matrix holding the _initial_ values of the 
       rate (i.e., inverse scale) parameters for the gamma prior 
       distributions over the v parameters.  i.e., beta[d,n] is the rate
       parameter governing v[d,n]. Introduced in eqn (16).
       NB:  this is the initial value beta, which is updated upon iteration.
       It is not (necessarily) the same as the hyperparameter beta0, which
       is unchanged by iteration.
  }
  \item{c}{
       a vector with D components holding the _initial_ values of the
       parameters of the Dirichlet distribution over the mixing 
       coefficients pi.  Introduced in eqn (19).
       NB: this is the initial value c, which is updated upon iteration.
       It is not (necessarily) the same as the hyperparameter c0, which 
       is unchanged by iteration.
  }
  \item{r}{
       the N x N.c matrix of initial responsibilities, with r[n, nc] giving the
       probability that item n belongs to component nc
  }
  \item{kmeans.clusters}{
       an N-vector giving the assignment of each of the N
       items to a cluster, as determined by kmeans.
  }
  \item{kmeans.centers}{
       an N.c x D matrix holding the centers of the N.c
       clusters/components determined by kmeans
  }
}

