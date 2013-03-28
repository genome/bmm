\name{bmm}
\alias{bmm}
\title{
  bmm:  Fits a Beta mixture with fixed number of components using a variational Bayesian approach
}
\description{
   Use a variational Bayesian approach to fit a mixture of Beta 
   distributions to proportion data, dropping any components/clusters that
   have small probability/mass.

   NB: Clustering is first run to convergence using bmm.fixed.num.components.
   If any components have probability less than pi.threshold, they are
   discarded and the clustering is run to convergence again.

   The core implements the derivations described in

   Bayesian Estimation of Beta Mixture Models with Variational
   Inference.  Ma and Leijon.  IEEE Transactions on Pattern Analysis
   and Machine Intelligence (2011) 33: 2160-2173.

   and

   Variational Learning for Finite Dirichlet Mixture Models and
   Applications.  Fan, Bouguila, and Ziou.  IEEE Transactions on
   Neural Networks and Learning Systems (2012) 23: 762-774.

   Notation and references here follow that used in Ma and Leijon.
}
\usage{
   bmm(X, N.c, r, mu, alpha, nu, beta, c, 
                            mu0, alpha0, nu0, beta0, c0, 
                            convergence.threshold = 10^-4,
                            max.iterations = 10000,
                            verbose = 0,
                            pi.convergence = 10^-2)
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
  \item{r}{
       the N x N.c matrix of initial responsibilities, with r[n, nc] giving the
       probability that item n belongs to component nc
  }
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
  \item{mu0, alpha0, nu0, beta0, c0}{  
       the hyperparameters corresponding to the
       above initial values (and with the same
       respective matrix/vector dimensionality).
  }
  \item{convergence.threshold}{
       minimum absolute difference between mixing
       coefficient (expected) values across consecutive
       iterations to reach converge.
  }
  \item{max.iterations}{
       maximum number of iterations to attempt
  }
  \item{verbose}{
       output progress in terms of mixing coefficient (expected) values if 1.
  }
  \item{pi.threshold}{
       discard any cluster with weight/mixing coefficient less
       than pi.threshold _following_ convergence.  
  }
}

\details{}

\value{
  A list with the following entries:
  \item{retVal}{
       0 indicates successful convergence; -1 indicates a failure
       to converge.
  }
  \item{mu}{
       a D x N.c matrix holding the _converged final_ values of the 
       shape parameters for the gamma prior distributions over the 
       u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
       Introduced in eqn (15).
  }
  \item{alpha}{
       a D x N.c matrix holding the _converged final_ values of the 
       rate (i.e., inverse scale) parameters for the gamma prior 
       distributions over the u parameters.  i.e., mu[d,n] is the rate
       parameter governing u[d,n]. Introduced in eqn (15).
  }
  \item{nu}{
       a D x N.c matrix holding the _converged final_ values of the 
       shape parameters for the gamma prior distributions over the 
       v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
       Introduced in eqn (16).
  }
  \item{beta}{
       a D x N.c matrix holding the _converged final_ values of the 
       rate (i.e., inverse scale) parameters for the gamma prior 
       distributions over the v parameters.  i.e., beta[d,n] is the rate
       parameter governing v[d,n]. Introduced in eqn (16).
  }
  \item{c}{
       a vector with D components holding the _converged final_ values of the
       parameters of the Dirichlet distribution over the mixing 
       coefficients pi.  Introduced in eqn (19).
  }
  \item{r}{
       the N x N.c matrix of responsibilities, with r[n, nc] giving the
       probability that item n belongs to component nc
  }
  \item{num.iterations}{
       the number of iterations required to reach convergence.
  }
  \item{ln.rho}{
       an N x N.c matrix holding the ln[rho], as defined in eqn (32).
  }
  \item{E.lnu}{
       the D x N.c matrix holding the values E_u[ln u], defined following
       eqn (51).
  }
  \item{E.lnv}{
       the D x N.c matrix holding the values E_v[ln v], defined following
       eqn (51).
  }
  \item{E.lnpi}{
       the D-vector holding the values E[ln pi], defined following 
       eqn (51).
  }
  \item{E.pi}{
       the D-vector holding the values E[pi], i.e., the expected values
       of the mixing coefficients, defined in eqn (53).
  }
  \item{E.quadratic.u}{
       the D x N.c matrix holding the values 
       E_u[(ln u - ln u^bar)^2] defined following eqn (51).
  }
  \item{E.quadratic.v}{
       the D x N.c matrix holding the values 
       E_v[(ln v - ln v^bar)^2] defined following eqn (51).
  }
  \item{ubar}{
       the D x N.c matrix holding values ubar = mu/alpha defined
       following eqn (51).
  }
  \item{vbar}{  
       the D x N.c matrix holding values vbar = nu/beta defined
       following eqn (51).
  }
}
\examples{}

