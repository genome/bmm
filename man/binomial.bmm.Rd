\name{binomial.bmm}
\alias{binomial.bmm}
\title{
  binomial.bmm: Fits a binomial mixture using a variational Bayesian approach
}
\description{
  Use a variational Bayesian approach to fit a mixture of binomial
  distributions to count data, dropping any components/clusters that
  have small probability/mass.
  
  NB: Clustering is first run to convergence using 
  binomial.bmm.fixed.num.components.
  If any components have probability less than pi.threshold, they are
  discarded and the clustering is run to convergence again.
}
\usage{
  binomial.bmm(successes, total.trials, N.c, r, a, b, alpha,
               a0, b0, alpha0, convergence.threshold = 10^-4, 
               max.iterations = 10000, verbose = 0, pi.threshold = 10^-2)
}
\arguments{
  \item{successes}{
    an N x D matrix with rows holding the number of successes
    in each dimension.
    All entries are assumed to be counts--not proportions.
  }
  \item{total.trials}{
    an N x D matrix with rows holding the number of trials
    in each dimension.
    All entries are assumed to be counts--not proportions.
  }
  \item{N.c}{
    the number of components/clusters to attempt
  }
  \item{r}{
    the N x N.c matrix of initial responsibilities, with r[n, nc] giving
    the probability that item n belongs to component nc
  }
  \item{a}{

  }
  \item{b}{

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
  \item{a0, b0, alpha0}{
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
    output progress in terms of mixing coefficient (expected) values if 1
  }
  \item{pi.threshold}{

  }
}

\value{
  A list with the following entries:
  \item{retVal}{
    0 indicates successful convergence; -1 indicates a failure to converge.
  }
  \item{a}{

  }
  \item{b}{

  }
  \item{alpha}{
    a D x N.c matrix holding the _converged final_ values of the 
    rate (i.e., inverse scale) parameters for the gamma prior 
    distributions over the u parameters.  i.e., mu[d,n] is the rate
    parameter governing u[d,n]. Introduced in eqn (15).
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
  \item{E.lnpi}{
    the D-vector holding the values E[ln pi], defined following eqn (51).
  }
  \item{E.pi}{
    the D-vector holding the values E[pi], i.e., the expected values
    of the mixing coefficients, defined in eqn (53).
  }
}
