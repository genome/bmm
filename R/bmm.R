suppressPackageStartupMessages(library(ggplot2))

####--------------------------------------------------------------
## bmm.fixed.num.components: Use a variational Bayesian approach to fit 
## a mixture of Beta distributions to proportion data, without dropping
## any components/clusters.  To instead automatically determine the number of
## components, use bmm, which invokes this function.
## This implements the derivations described in
##
## Bayesian Estimation of Beta Mixture Models with Variational
## Inference.  Ma and Leijon.  IEEE Transactions on Pattern Analysis
## and Machine Intelligence (2011) 33: 2160-2173.
##
## and
##
## Variational Learning for Finite Dirichlet Mixture Models and
## Applications.  Fan, Bouguila, and Ziou.  IEEE Transactions on
## Neural Networks and Learning Systems (2012) 23: 762-774.
##
## Notation and references here follow that used in Ma and Leijon.
##
## Inputs:
## X:  an N x D matrix with rows being the items to cluster.
##     All entries are assumed to be proportions (i.e., between 0 and 1).
##     Notice that there are no summation restrictions--i.e., proportions
##     do not sum to unity across an item's dimensions.
## N.c:  the number of components/clusters to attempt
## r:  the N x N.c matrix of initial responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## mu:  a D x N.c matrix holding the _initial_ values of the 
##      shape parameters for the gamma prior distributions over the 
##      u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
##      NB:  this is the initial value mu, which is updated upon iteration.
##      It is not (necessarily) the same as the hyperparameter mu0, which
##      is unchanged by iteration.
##      Introduced in eqn (15).
## alpha:  a D x N.c matrix holding the _initial_ values of the 
##         rate (i.e., inverse scale) parameters for the gamma prior 
##         distributions over the u parameters.  i.e., mu[d,n] is the rate
##         parameter governing u[d,n]. Introduced in eqn (15).
##         NB:  this is the initial value alpha, which is updated upon iteration.
##         It is not (necessarily) the same as the hyperparameter alpha0, which
##         is unchanged by iteration.
## nu:  a D x N.c matrix holding the _initial_ values of the 
##      shape parameters for the gamma prior distributions over the 
##      v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
##      Introduced in eqn (16).
##      NB:  this is the initial value nu, which is updated upon iteration.
##      It is not (necessarily) the same as the hyperparameter nu0, which
##      is unchanged by iteration.
## beta:  a D x N.c matrix holding the _initial_ values of the 
##        rate (i.e., inverse scale) parameters for the gamma prior 
##        distributions over the v parameters.  i.e., beta[d,n] is the rate
##        parameter governing v[d,n]. Introduced in eqn (16).
##        NB:  this is the initial value beta, which is updated upon iteration.
##        It is not (necessarily) the same as the hyperparameter beta0, which
##        is unchanged by iteration.
## c:  a vector with D components holding the _initial_ values of the
##     parameters of the Dirichlet distribution over the mixing 
##     coefficients pi.  Introduced in eqn (19).
##     NB: this is the initial value c, which is updated upon iteration.
##     It is not (necessarily) the same as the hyperparameter c0, which 
##     is unchanged by iteration.
## mu0, alpha0, nu0, beta0, c0:  the hyperparameters corresponding to the
##                               above initial values (and with the same
##                               respective matrix/vector dimensionality).
## convergence.threshold:  minimum absolute difference between mixing
##                         coefficient (expected) values across consecutive
##                         iterations to reach converge.
## max.iterations:  maximum number of iterations to attempt
## verbose:  output progress in terms of mixing coefficient (expected) values
##           if 1.
## Outputs: a list with the following entries
## retVal:  0 indicates successful convergence; -1 indicates a failure
##          to converge.
## mu:  a D x N.c matrix holding the _converged final_ values of the 
##      shape parameters for the gamma prior distributions over the 
##      u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
##      Introduced in eqn (15).
## alpha:  a D x N.c matrix holding the _converged final_ values of the 
##         rate (i.e., inverse scale) parameters for the gamma prior 
##         distributions over the u parameters.  i.e., mu[d,n] is the rate
##         parameter governing u[d,n]. Introduced in eqn (15).
## nu:  a D x N.c matrix holding the _converged final_ values of the 
##      shape parameters for the gamma prior distributions over the 
##      v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
##      Introduced in eqn (16).
## beta:  a D x N.c matrix holding the _converged final_ values of the 
##        rate (i.e., inverse scale) parameters for the gamma prior 
##        distributions over the v parameters.  i.e., beta[d,n] is the rate
##        parameter governing v[d,n]. Introduced in eqn (16).
## c:  a vector with D components holding the _converged final_ values of the
##     parameters of the Dirichlet distribution over the mixing 
##     coefficients pi.  Introduced in eqn (19).
## r:  the N x N.c matrix of responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## num.iterations:  the number of iterations required to reach convergence.
## ln.rho:  an N x N.c matrix holding the ln[rho], as defined in eqn (32).
## E.lnu:  the D x N.c matrix holding the values E_u[ln u], defined following
##         eqn (51).
## E.lnv:  the D x N.c matrix holding the values E_v[ln v], defined following
##         eqn (51).
## E.lnpi:  the D-vector holding the values E[ln pi], defined following 
##          eqn (51).
## E.pi:  the D-vector holding the values E[pi], i.e., the expected values
##        of the mixing coefficients, defined in eqn (53).
## E.quadratic.u:  the D x N.c matrix holding the values 
##                 E_u[(ln u - ln u^bar)^2] defined following eqn (51).
## E.quadratic.v:  the D x N.c matrix holding the values 
##                 E_v[(ln v - ln v^bar)^2] defined following eqn (51).
## ubar:  the D x N.c matrix holding values ubar = mu/alpha defined
##        following eqn (51).
## vbar:  the D x N.c matrix holding values vbar = nu/beta defined
##        following eqn (51).
bmm.fixed.num.components <- function(X, N.c, r, mu, alpha, nu, beta, c, E.pi, mu0, alpha0, nu0, beta0, c0, convergence.threshold = 10^-4, max.iterations = 10000, verbose = 0)
{
  ubar <- mu / alpha
  vbar <- nu / beta
  
  N <- dim(X)[1]
  D <- dim(X)[2]

  E.pi.prev <- rep(0, N.c)
  E.pi.prev <- E.pi
  
  iteration <- 0

  # Apply variational bayesian approach to Beta mixture modeling
  # until convergence
  while(TRUE) {  

    #expected.means <- ubar / ( ubar + vbar )
    #cat("Expected means:\n")
    #print(expected.means)
    
    # E_u[ln u] defined following eqn (51).
    E.lnu <- digamma(mu) - log(alpha)

    # E_v[ln v] defined following eqn (51).
    E.lnv <- digamma(nu) - log(beta)
    
    # E[ln pi_i] defined following eqn (51).
    E.lnpi <- digamma(c) - digamma(sum(c))
  
    # E[pi_i] defined in eqn (53).
    #E.pi <- ( c0 + colSums(r, na.rm=TRUE) ) / ( sum(c0) + N )

    # E_u[(ln u - ln u^bar)^2] defined following eqn (51).
    E.quadratic.u <- ( ( digamma(mu) - log(mu) )^2 ) + trigamma(mu)
  
    # E_v[(ln v - ln v^bar)^2] defined following eqn (51).
    E.quadratic.v <- ( ( digamma(nu) - log(nu) )^2 ) + trigamma(nu)
    
    # Define some temporary values to be used below.
    lg.u.v <- lgamma(ubar + vbar)
    lg.u <- lgamma(ubar)
    lg.v <- lgamma(vbar)
    dig.u.v <- digamma(ubar + vbar)
    dig.u <- digamma(ubar)
    dig.v <- digamma(vbar)
    trig.u.v <- trigamma(ubar + vbar)
    trig.u <- trigamma(ubar)
    trig.v <- trigamma(vbar)
    E.lnv.logvbar <- E.lnv - log(vbar)
    E.lnu.logubar <- E.lnu - log(ubar)
    u.logx <- log(X) %*% (ubar-1)
    v.logx <- log(1-X) %*% (vbar-1)

    # ln[rho] defined in eqn (32).
    tmp <- colSums((lg.u.v - lg.u - lg.v) + (ubar * (dig.u.v - dig.u) * E.lnu.logubar) + (vbar * (dig.u.v - dig.v) * E.lnv.logvbar) + 0.5 * (ubar^2 * (trig.u.v - trig.u) * E.quadratic.u) + 0.5 * (vbar^2 * (trig.u.v - trig.v) * E.quadratic.v) + ubar * vbar * trig.u.v * E.lnu.logubar * E.lnv.logvbar)
    tmp.matrix <- matrix(data=tmp, nrow=N, ncol=N.c, byrow=TRUE)
    ln.rho <- matrix(data=E.lnpi, nrow=N, ncol=N.c, byrow=TRUE)

    ln.rho <- ln.rho + u.logx + v.logx + tmp.matrix

    # Responsibilities r defined in eqn (31).
    r <- matrix(data = 0, nrow=N, ncol=N.c)
    for(n in 1:N) {
      if(any(is.na(ln.rho[n,]))) {
        r[n,] <- rep(NA, N.c)
        next
      }
      row.sum <- log(sum(exp(ln.rho[n,] - max(ln.rho[n,])))) + max(ln.rho[n,], na.rm=TRUE)
      for(k in 1:N.c) { r[n,k] = exp(ln.rho[n,k] - row.sum) }
    }

    r <- r + 10^-9

    # r and X will have NAs for items that have been removed.
    # Substitute 0's for these NAs.
    r.na.rows <- is.na(rowSums(r))
    X.na.rows <- is.na(rowSums(X))
    if(any(is.na(r.na.rows))) {
      print(r.na.rows)
      print(X.na.rows)
    }
    if(any(r.na.rows != X.na.rows)) {
      print(r)
      print(X)
      stop(sprintf("NA's inconsistent between r and X matrices\n"))
    }
    r.no.na <- r[!r.na.rows,]
    X.no.na <- X[!X.na.rows,]

    # E[pi_i] defined in eqn (53).
    E.pi <- ( c0 + colSums(r, na.rm=TRUE) ) / ( sum(c0) + N )
    
    # Update alpha as defined in eqn (49).
    alpha <- alpha0 - t(t(r.no.na) %*% log(X.no.na))

    # Update beta as defined in eqn (51).
    # beta <- beta0 - t(t(r.no.na) %*% log(1-X.no.na))
    beta <- beta0 - t(t(r.no.na) %*% log(1-X.no.na))
  
    # Update mu as defined in eqn (48).
    v.trig.E.log <- vbar * trig.u.v * E.lnv.logvbar
    r.colsums <- matrix(data=colSums(r, na.rm=TRUE), nrow=D, ncol=N.c, byrow=TRUE)
    mu <- mu0 + r.colsums * ubar * ( dig.u.v - dig.u + v.trig.E.log )
     
    # Update nu as defined in eqn (50).
    u.trig.E.log <- ubar * trig.u.v * E.lnu.logubar
    nu <- nu0 + r.colsums * vbar * ( dig.u.v - dig.v + u.trig.E.log )

    # Update c as defined in eqn (47).
    c <- c0 + colSums(r, na.rm=TRUE)
  
    # E-step (equations following eqn 51) on page 8
    ubar <- mu / alpha
    vbar <- nu / beta
    
    if(verbose > 0) {
      print(E.pi)
    } 

    # Convergence test
    if (( all(abs(E.pi - E.pi.prev) < convergence.threshold) )) {
      break
    }
        
    E.pi.prev <- E.pi
    
    iteration <- iteration + 1

    if(iteration >= max.iterations) {
      cat("Exceeded ", max.iterations, " iterations\n")
      retList <- list("retVal" = -1, "mu" = mu, "alpha" = alpha, "nu" = nu, "beta" = beta, "c" = c, "r" = r, "num.iterations" = iteration, "ln.rho" = ln.rho, "E.lnu" = E.lnu, "E.lnv" = E.lnv, "E.lnpi" = E.lnpi, "E.pi" = E.pi, "E.quadratic.u" = E.quadratic.u, "E.quadratic.v" = E.quadratic.v, "ubar" = ubar, "vbar" = vbar)      
      return(retList)
    }
    
  } # End inner while(TRUE)

  retList <- list("retVal" = 0, "mu" = mu, "alpha" = alpha, "nu" = nu, "beta" = beta, "c" = c, "r" = r, "num.iterations" = iteration, "ln.rho" = ln.rho, "E.lnu" = E.lnu, "E.lnv" = E.lnv, "E.lnpi" = E.lnpi, "E.pi" = E.pi, "E.quadratic.u" = E.quadratic.u, "E.quadratic.v" = E.quadratic.v, "ubar" = ubar, "vbar" = vbar)

  return(retList)
} # End bmm.fixed.num.components function

####--------------------------------------------------------------
## bmm: Use a variational Bayesian approach to fit a mixture of Beta 
## distributions to proportion data, dropping any components/clusters that
## have small probability/mass.
##
## NB: Clustering is first run to convergence using bmm.fixed.num.components.
## If any components have probability less than pi.threshold, they are
## discarded and the clustering is run to convergence again.
##
## The core implements the derivations described in
##
## Bayesian Estimation of Beta Mixture Models with Variational
## Inference.  Ma and Leijon.  IEEE Transactions on Pattern Analysis
## and Machine Intelligence (2011) 33: 2160-2173.
##
## and
##
## Variational Learning for Finite Dirichlet Mixture Models and
## Applications.  Fan, Bouguila, and Ziou.  IEEE Transactions on
## Neural Networks and Learning Systems (2012) 23: 762-774.
##
## Notation and references here follow that used in Ma and Leijon.
##
## Inputs:
## X:  an N x D matrix with rows being the items to cluster.
##     All entries are assumed to be proportions (i.e., between 0 and 1).
##     Notice that there are no summation restrictions--i.e., proportions
##     do not sum to unity across an item's dimensions.
## N.c:  the number of components/clusters to attempt
## r:  the N x N.c matrix of initial responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## mu:  a D x N.c matrix holding the _initial_ values of the 
##      shape parameters for the gamma prior distributions over the 
##      u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
##      NB:  this is the initial value mu, which is updated upon iteration.
##      It is not (necessarily) the same as the hyperparameter mu0, which
##      is unchanged by iteration.
##      Introduced in eqn (15).
## alpha:  a D x N.c matrix holding the _initial_ values of the 
##         rate (i.e., inverse scale) parameters for the gamma prior 
##         distributions over the u parameters.  i.e., mu[d,n] is the rate
##         parameter governing u[d,n]. Introduced in eqn (15).
##         NB:  this is the initial value alpha, which is updated upon iteration.
##         It is not (necessarily) the same as the hyperparameter alpha0, which
##         is unchanged by iteration.
## nu:  a D x N.c matrix holding the _initial_ values of the 
##      shape parameters for the gamma prior distributions over the 
##      v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
##      Introduced in eqn (16).
##      NB:  this is the initial value nu, which is updated upon iteration.
##      It is not (necessarily) the same as the hyperparameter nu0, which
##      is unchanged by iteration.
## beta:  a D x N.c matrix holding the _initial_ values of the 
##        rate (i.e., inverse scale) parameters for the gamma prior 
##        distributions over the v parameters.  i.e., beta[d,n] is the rate
##        parameter governing v[d,n]. Introduced in eqn (16).
##        NB:  this is the initial value beta, which is updated upon iteration.
##        It is not (necessarily) the same as the hyperparameter beta0, which
##        is unchanged by iteration.
## c:  a vector with D components holding the _initial_ values of the
##     parameters of the Dirichlet distribution over the mixing 
##     coefficients pi.  Introduced in eqn (19).
##     NB: this is the initial value c, which is updated upon iteration.
##     It is not (necessarily) the same as the hyperparameter c0, which 
##     is unchanged by iteration.
## mu0, alpha0, nu0, beta0, c0:  the hyperparameters corresponding to the
##                               above initial values (and with the same
##                               respective matrix/vector dimensionality).
## convergence.threshold:  minimum absolute difference between mixing
##                         coefficient (expected) values across consecutive
##                         iterations to reach converge.
## max.iterations:  maximum number of iterations to attempt
## verbose:  output progress in terms of mixing coefficient (expected) values
##           if 1.
## pi.threshold:  discard any cluster with weight/mixing coefficient less
##                than pi.threshold _following_ convergence.  
## Outputs: a list with the following entries
## retVal:  0 indicates successful convergence; -1 indicates a failure
##          to converge.
## mu:  a D x N.c matrix holding the _converged final_ values of the 
##      shape parameters for the gamma prior distributions over the 
##      u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
##      Introduced in eqn (15).
## alpha:  a D x N.c matrix holding the _converged final_ values of the 
##         rate (i.e., inverse scale) parameters for the gamma prior 
##         distributions over the u parameters.  i.e., mu[d,n] is the rate
##         parameter governing u[d,n]. Introduced in eqn (15).
## nu:  a D x N.c matrix holding the _converged final_ values of the 
##      shape parameters for the gamma prior distributions over the 
##      v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
##      Introduced in eqn (16).
## beta:  a D x N.c matrix holding the _converged final_ values of the 
##        rate (i.e., inverse scale) parameters for the gamma prior 
##        distributions over the v parameters.  i.e., beta[d,n] is the rate
##        parameter governing v[d,n]. Introduced in eqn (16).
## c:  a vector with D components holding the _converged final_ values of the
##     parameters of the Dirichlet distribution over the mixing 
##     coefficients pi.  Introduced in eqn (19).
## r:  the N x N.c matrix of responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## num.iterations:  the number of iterations required to reach convergence.
## ln.rho:  an N x N.c matrix holding the ln[rho], as defined in eqn (32).
## E.lnu:  the D x N.c matrix holding the values E_u[ln u], defined following
##         eqn (51).
## E.lnv:  the D x N.c matrix holding the values E_v[ln v], defined following
##         eqn (51).
## E.lnpi:  the D-vector holding the values E[ln pi], defined following 
##          eqn (51).
## E.pi:  the D-vector holding the values E[pi], i.e., the expected values
##        of the mixing coefficients, defined in eqn (53).
## E.quadratic.u:  the D x N.c matrix holding the values 
##                 E_u[(ln u - ln u^bar)^2] defined following eqn (51).
## E.quadratic.v:  the D x N.c matrix holding the values 
##                 E_v[(ln v - ln v^bar)^2] defined following eqn (51).
## ubar:  the D x N.c matrix holding values ubar = mu/alpha defined
##        following eqn (51).
## vbar:  the D x N.c matrix holding values vbar = nu/beta defined
##        following eqn (51).

bmm <- function(X, N.c, r, mu, alpha, nu, beta, c, mu0, alpha0, nu0, beta0, c0, convergence.threshold = 10^-4, max.iterations = 10000, verbose = 0, pi.threshold = 10^-2)
{

  total.iterations <- 0 
  D <- dim(X)[2]

  E.pi <- rep(1/N.c, N.c)
  
  while(TRUE) {

    if(verbose){
      print(r)
    }
    
    bmm.res <- bmm.fixed.num.components(X, N.c, r, mu, alpha, nu, beta, c, E.pi, mu0, alpha0, nu0, beta0, c0, convergence.threshold = 10^-4, max.iterations = 10000, verbose = verbose)
    if(bmm.res$retVal != 0) {
      stop("Failed to converge!\n")
    }

    mu <- bmm.res$mu
    alpha <- bmm.res$alpha
    nu <- bmm.res$nu
    beta <- bmm.res$beta
    c <- bmm.res$c
    E.pi <- bmm.res$E.pi

    total.iterations <- total.iterations + bmm.res$num.iterations

    ln.rho <- bmm.res$ln.rho
    E.lnu <- bmm.res$E.lnu
    E.lnv <- bmm.res$E.lnv
    E.lnpi <- bmm.res$E.lnpi
    E.quadratic.u <- bmm.res$E.quadratic.u
    E.quadratic.v <- bmm.res$E.quadratic.v
    ubar <- bmm.res$ubar
    vbar <- bmm.res$vbar

    do.inner.iteration <- FALSE

    apply.min.items.condition <- TRUE

    if((apply.min.items.condition == TRUE) & (N.c > 1)) {
      non.zero.indices <- E.pi > pi.threshold
      N = length(X[,1])
      clusters <- rep(0, N)
      for(n in 1:N) {
        max.cluster <- 0
        max.assignment <- -1
        for(k in 1:N.c) {
          if ( r[n,k] > max.assignment ) {
            max.assignment <- r[n,k]
            max.cluster <- k
          }
        }
        clusters[n] <- max.cluster
      }

      if ( any(non.zero.indices==FALSE) ) {

        do.inner.iteration <- TRUE

        numeric.indices <- (1:N.c)

        E.pi <- E.pi[non.zero.indices]
        E.lnpi <- E.lnpi[non.zero.indices]

        N.c <- length(E.pi)
        c <- c[non.zero.indices]
        c0 <- c0[non.zero.indices]
        r <- matrix(r[,non.zero.indices], nrow=N, ncol=N.c)
        ln.rho <- matrix(ln.rho[,non.zero.indices], nrow=N, ncol=N.c)
        # Need to renormalize r--do it gently.
        for(n in 1:N) {
           if(any(is.na(ln.rho[n,]))) {
             r[n,] <- rep(NA, N.c)
             next
           }
          
           row.sum <- log(sum(exp(ln.rho[n,] - max(ln.rho[n,])))) + max(ln.rho[n,])
           for(k in 1:N.c) { r[n,k] = exp(ln.rho[n,k] - row.sum) }
        }
        mu <- matrix(mu[,non.zero.indices], nrow=D, ncol=N.c)
        nu <- matrix(nu[,non.zero.indices], nrow=D, ncol=N.c)
        mu0 <- matrix(mu0[,non.zero.indices], nrow=D, ncol=N.c)
        nu0 <- matrix(nu0[,non.zero.indices], nrow=D, ncol=N.c)
        alpha <- matrix(alpha[,non.zero.indices], nrow=D, ncol=N.c)
        beta <- matrix(beta[,non.zero.indices], nrow=D, ncol=N.c)
        alpha0 <- matrix(alpha0[,non.zero.indices], nrow=D, ncol=N.c)
        beta0 <- matrix(beta0[,non.zero.indices], nrow=D, ncol=N.c)
      }
    } # End apply.min.items.condition

    if(do.inner.iteration == FALSE) { break }

  }

  retList <- list("retVal" = 0, "mu" = mu, "alpha" = alpha, "nu" = nu, "beta" = beta, "c" = c, "r" = r, "num.iterations" = total.iterations, "ln.rho" = ln.rho, "E.lnu" = E.lnu, "E.lnv" = E.lnv, "E.lnpi" = E.lnpi, "E.pi" = E.pi, "E.quadratic.u" = E.quadratic.u, "E.quadratic.v" = E.quadratic.v, "ubar" = ubar, "vbar" = vbar)

  return(retList)

} # End bmm function

####--------------------------------------------------------------
## init.bmm.hyperparameters:  Initialize the bmm hyperparameters (to be
##                            passed to bmm)
##
## NB:  This provides what should be a generally reasonable initialization of
##      hyperparameters.  However, better results may be obtained by
##      tuning these in an application-specific manner.
##
## Inputs:
## X:  an N x D matrix with rows being the items to cluster.
##     All entries are assumed to be proportions (i.e., between 0 and 1).
##     Notice that there are no summation restrictions--i.e., proportions
##     do not sum to unity across an item's dimensions.
## N.c:  the number of components/clusters to attempt
## Outputs:
## mu0:  a D x N.c matrix holding the hyperparameter values of the 
##       shape parameters for the gamma prior distributions over the 
##       u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
##       Introduced in eqn (15).
## alpha0:  a D x N.c matrix holding the hyperparameter values of the 
##          rate (i.e., inverse scale) parameters for the gamma prior 
##          distributions over the u parameters.  i.e., mu[d,n] is the rate
##          parameter governing u[d,n]. Introduced in eqn (15).
## nu0:  a D x N.c matrix holding the hyperparameter values of the 
##       shape parameters for the gamma prior distributions over the 
##       v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
##       Introduced in eqn (16).
## beta0:  a D x N.c matrix holding the hyperparameter values of the 
##         rate (i.e., inverse scale) parameters for the gamma prior 
##         distributions over the v parameters.  i.e., beta[d,n] is the rate
##         parameter governing v[d,n]. Introduced in eqn (16).
## c0:  a vector with D components holding the hyperparameter values of the
##      parameters of the Dirichlet distribution over the mixing 
##      coefficients pi.  Introduced in eqn (19).
init.bmm.hyperparameters <- function(X, N.c) 
{
  N <- dim(X)[1]
  D <- dim(X)[2]

  # I. Choose the initial parameteres C_10, ..., C_I0 for the
  # Dirichlet distribution.
  # Follow Ma and Leijon in setting c_i0 = 0.001 for all i,
  # which ensures that the (numbers of) components are determined primarily
  # from the data, not from the prior.
  # c is alpha in Fan's notation.
  c0 <- rep(0.001, N.c)

  # II.  Choose the initial parameters (element-wise) alpha0 > 0,
  # beta0 > 0, mu0 > 0.6156, nu0 > 0.6156.  Furthermore,
  # mu0 / alpha0 and nu0 / beta0 should be greater than 1.
  
  # Here we choose alpha0 and beta0 the same ...
  # By choosing the alpha parameter small, we make the prior flat(ter)
  # alpha0/beta0 is v in Fan's notation
  alpha0 <- matrix(data=0.005, nrow=D, ncol=N.c)
  #alpha0 <- matrix(data=0.0075, nrow=D, ncol=N.c)    
  #alpha0 <- matrix(data=0.01, nrow=D, ncol=N.c)  
  beta0 <- alpha0

  # ... and mu0 and nu0 the same.
  # NB:  with choice mu=1, the gamma prior degenerates to an
  # exponential distribution.
  # mu/nu is u in Fan's notation
  mu0 <- matrix(data=1, nrow=D, ncol=N.c)
  nu0 <- mu0

  retList <- list("mu0" = mu0, "alpha0" = alpha0, "nu0" = nu0, "beta0" = beta0, "c0" = c0)
  return(retList)

} # End init.bmm.hyperparameters


####--------------------------------------------------------------
## init.bmm.parameters:  Initialize the bmm parameters (to be
##                       passed to bmm)
##
## Initialize parameters such that expected proportions have the values
## determined by an initial k-means clustering.
##
## NB:  This provides what should be a generally reasonable initialization of
##      hyperparameters.  However, better results may be obtained by
##      tuning these in an application-specific manner.
##
##
## Inputs:
## X:  an N x D matrix with rows being the items to cluster.
##     All entries are assumed to be proportions (i.e., between 0 and 1).
##     Notice that there are no summation restrictions--i.e., proportions
##     do not sum to unity across an item's dimensions.
## N.c:  the number of components/clusters to attempt
## mu0:  a D x N.c matrix holding the hyperparameter values of the 
##       shape parameters for the gamma prior distributions over the 
##       u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
##       Introduced in eqn (15).
## alpha0:  a D x N.c matrix holding the hyperparameter values of the 
##          rate (i.e., inverse scale) parameters for the gamma prior 
##          distributions over the u parameters.  i.e., mu[d,n] is the rate
##          parameter governing u[d,n]. Introduced in eqn (15).
## nu0:  a D x N.c matrix holding the hyperparameter values of the 
##       shape parameters for the gamma prior distributions over the 
##       v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
##       Introduced in eqn (16).
## beta0:  a D x N.c matrix holding the hyperparameter values of the 
##         rate (i.e., inverse scale) parameters for the gamma prior 
##         distributions over the v parameters.  i.e., beta[d,n] is the rate
##         parameter governing v[d,n]. Introduced in eqn (16).
## c0:  a vector with D components holding the hyperparameter values of the
##      parameters of the Dirichlet distribution over the mixing 
##      coefficients pi.  Introduced in eqn (19).
##
## Outputs:
## mu:  a D x N.c matrix holding the _initial_ values of the 
##      shape parameters for the gamma prior distributions over the 
##      u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
##      NB:  this is the initial value mu, which is updated upon iteration.
##      It is not (necessarily) the same as the hyperparameter mu0, which
##      is unchanged by iteration.
##      Introduced in eqn (15).
## alpha:  a D x N.c matrix holding the _initial_ values of the 
##         rate (i.e., inverse scale) parameters for the gamma prior 
##         distributions over the u parameters.  i.e., mu[d,n] is the rate
##         parameter governing u[d,n]. Introduced in eqn (15).
##         NB:  this is the initial value alpha, which is updated upon iteration.
##         It is not (necessarily) the same as the hyperparameter alpha0, which
##         is unchanged by iteration.
## nu:  a D x N.c matrix holding the _initial_ values of the 
##      shape parameters for the gamma prior distributions over the 
##      v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
##      Introduced in eqn (16).
##      NB:  this is the initial value nu, which is updated upon iteration.
##      It is not (necessarily) the same as the hyperparameter nu0, which
##      is unchanged by iteration.
## beta:  a D x N.c matrix holding the _initial_ values of the 
##        rate (i.e., inverse scale) parameters for the gamma prior 
##        distributions over the v parameters.  i.e., beta[d,n] is the rate
##        parameter governing v[d,n]. Introduced in eqn (16).
##        NB:  this is the initial value beta, which is updated upon iteration.
##        It is not (necessarily) the same as the hyperparameter beta0, which
##        is unchanged by iteration.
## c:  a vector with D components holding the _initial_ values of the
##     parameters of the Dirichlet distribution over the mixing 
##     coefficients pi.  Introduced in eqn (19).
##     NB: this is the initial value c, which is updated upon iteration.
##     It is not (necessarily) the same as the hyperparameter c0, which 
##     is unchanged by iteration.
## r:  the N x N.c matrix of initial responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## kmeans.clusters:  an N-vector giving the assignment of each of the N
##                   items to a cluster, as determined by kmeans.
## kmeans.centers:  an N.c x D matrix holding the centers of the N.c
##                  clusters/components determined by kmeans
init.bmm.parameters <- function(X, N.c, mu0, alpha0, nu0, beta0, c0, verbose=0)
{
  N <- dim(X)[1]
  D <- dim(X)[2]

  # Initialize the responsibilities r_ni with K-means.
  kmeans.out <- kmeans(X, N.c, nstart=1000)
  kmeans.clusters <- kmeans.out$cluster
  kmeans.centers <- kmeans.out$centers

  # Sort the centers to ease comparison across runs.  No we sort the
  # clustered results.
  # Besides, be careful!  This upset the relationship between kmeans.clusters
  # and kmeans.centers!
  #kmeans.centers <- kmeans.centers[do.call(order, lapply(1:ncol(kmeans.centers), function(i) kmeans.centers[, i])), ]
  if(verbose){
    cat("kmeans initialization:\n")
    write.table(kmeans.centers, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
  }
  # This logic (should) follow that of Fan in DMM.m in initializing
  # the model.
  # (1) Initialize responsibilities with k-means
  # (2) Set the parameters to the hyperparameters
  # (3) E step using those parameters (without updating responsibilities)
  # (4) M step to update parameters

  r <- matrix(data=0, nrow=N, ncol=N.c)
  for(i in 1:N) {
    r[i,kmeans.clusters[i]]<- 1
  }

  # Initialize parameters to hyperparameters
  mu <- mu0
  nu <- nu0
  beta <- beta0
  alpha <- alpha0  

  # Set initial alpha and beta (NB: not the prior hyperparameters) according
  # to the means found by kmeans.  Keep mu and nu at their
  # prior values.  To do so:
  # E[proportion] = ubar / ( ubar + vbar ),
  #   where ubar = mu / alpha and vbar = nu / beta
  # For now, keep beta at prior as well.
  # Uncommenting this also seems to work fairly well.
  alpha <- ( mu * beta * ( 1 - t(kmeans.centers) ) ) / ( t(kmeans.centers) * nu )

  c <- c0 + colSums(r, na.rm=TRUE)
  
  # E.pi <- ( c0 + colSums(r, na.rm=TRUE) ) / ( sum(c0) + N )
  # NB: E.pi is only used for iteration, to test for convergence.
  # When it is needed for a calculation, it is computed from the parameters.
  E.pi <- rep(0, N.c)  

  ubar <- mu / alpha
  vbar <- nu / beta

  # E_u[ln u] defined following eqn (51).
  E.lnu <- digamma(mu) - log(alpha)

  # E_v[ln v] defined following eqn (51).
  E.lnv <- digamma(nu) - log(beta)
    
  # Update c as defined in eqn (47).
  c <- c0 + colSums(r, na.rm=TRUE)

  # E[ln pi_i] defined following eqn (51).
  E.lnpi <- digamma(c) - digamma(sum(c))
  
  # E[pi_i] defined in eqn (53).
  #E.pi <- ( c0 + colSums(r, na.rm=TRUE) ) / ( sum(c0) + N )
  # NB: E.pi is only used for iteration, to test for convergence.
  # When it is needed for a calculation, it is computed from the parameters.
  E.pi <- rep(0, N.c)  
  
  # E_u[(ln u - ln u^bar)^2] defined following eqn (51).
  E.quadratic.u <- ( ( digamma(mu) - log(mu) )^2 ) + trigamma(mu)
  
  # E_v[(ln v - ln v^bar)^2] defined following eqn (51).
  E.quadratic.v <- ( ( digamma(nu) - log(nu) )^2 ) + trigamma(nu)
    
  # Define some temporary values to be used below.
  lg.u.v <- lgamma(ubar + vbar)
  lg.u <- lgamma(ubar)
  lg.v <- lgamma(vbar)
  dig.u.v <- digamma(ubar + vbar)
  dig.u <- digamma(ubar)
  dig.v <- digamma(vbar)
  trig.u.v <- trigamma(ubar + vbar)
  trig.u <- trigamma(ubar)
  trig.v <- trigamma(vbar)
  E.lnv.logvbar <- E.lnv - log(vbar)
  E.lnu.logubar <- E.lnu - log(ubar)
  u.logx <- log(X) %*% (ubar-1)
  v.logx <- log(1-X) %*% (vbar-1)

  # E[pi_i] defined in eqn (53).
  #E.pi <- ( c0 + colSums(r, na.rm=TRUE) ) / ( sum(c0) + N )
    
  # Update alpha as defined in eqn (49).
  alpha <- alpha0 - t(t(r) %*% log(X))

  # Update beta as defined in eqn (51).
  # beta <- beta0 - t(t(r) %*% log(1-X))
  beta <- beta0 - t(t(r) %*% log(1-X))
  
  # Update mu as defined in eqn (48).
  v.trig.E.log <- vbar * trig.u.v * E.lnv.logvbar
  r.colsums <- matrix(data=colSums(r, na.rm=TRUE), nrow=D, ncol=N.c, byrow=TRUE)
  mu <- mu0 + r.colsums * ubar * ( dig.u.v - dig.u + v.trig.E.log )
     
  # Update nu as defined in eqn (50).
  u.trig.E.log <- ubar * trig.u.v * E.lnu.logubar
  nu <- nu0 + r.colsums * vbar * ( dig.u.v - dig.v + u.trig.E.log )

  retList <- list("mu" = mu, "alpha" = alpha, "nu" = nu, "beta" = beta, "c" = c, "E.pi" = E.pi, "r" = r, "kmeans.centers" = kmeans.centers, "kmeans.clusters" = kmeans.clusters)
  return(retList)
  
  # Old logic below.
  
  # Set initial alpha and beta (NB: not the prior hyperparameters) according
  # to the means found by kmeans.  Keep mu and nu at their
  # prior values.  To do so:
  # E[proportion] = ubar / ( ubar + vbar ),
  #   where ubar = mu / alpha and vbar = nu / beta
  # For now, keep beta at prior as well.
  mu <- mu0
  nu <- nu0
  beta <- beta0
  alpha <- ( mu * beta * ( 1 - t(kmeans.centers) ) ) / ( t(kmeans.centers) * nu )
  
  c <- c0 + colSums(r, na.rm=TRUE)
  
  # E.pi <- ( c0 + colSums(r, na.rm=TRUE) ) / ( sum(c0) + N )
  # NB: E.pi is only used for iteration, to test for convergence.
  # When it is needed for a calculation, it is computed from the parameters.
  E.pi <- rep(0, N.c)  

  #retList <- list("mu" = mu, "alpha" = alpha, "nu" = nu, "beta" = beta, "c" = c, "E.pi" = E.pi, "r" = r, "kmeans.centers" = kmeans.centers, "kmeans.clusters" = kmeans.clusters)
  #return(retList)

  ubar <- mu / alpha
  vbar <- nu / beta

  E.lnu <- digamma(mu) - log(alpha)
  E.lnv <- digamma(nu) - log(beta)

  # Update alpha as defined in eqn (49).
  alpha <- alpha0 - t(t(r) %*% log(X))

  # Update beta as defined in eqn (51).
  beta <- beta0 - t(t(r) %*% log(1-X))
  
  ubar <- mu / alpha
  vbar <- nu / beta

  E.lnu <- digamma(mu) - log(alpha)
  E.lnv <- digamma(nu) - log(beta)


  trig.u.v <- trigamma(ubar + vbar)
  E.lnv.logvbar <- E.lnv - log(vbar)
  E.lnu.logubar <- E.lnu - log(ubar)
  dig.u.v <- digamma(ubar + vbar)
  dig.u <- digamma(ubar)
  dig.v <- digamma(vbar)

  # Update mu as defined in eqn (48).
  v.trig.E.log <- vbar * trig.u.v * E.lnv.logvbar
  r.colsums <- matrix(data=colSums(r, na.rm=TRUE), nrow=D, ncol=N.c, byrow=TRUE)
  mu <- mu0 + r.colsums * ubar * ( dig.u.v - dig.u + v.trig.E.log )
     
  # Update nu as defined in eqn (50).
  u.trig.E.log <- ubar * trig.u.v * E.lnu.logubar
  nu <- nu0 + r.colsums * vbar * ( dig.u.v - dig.v + u.trig.E.log )

  c <- c0 + colSums(r, na.rm=TRUE) 
  
  retList <- list("mu" = mu, "alpha" = alpha, "nu" = nu, "beta" = beta, "c" = c, "E.pi" = E.pi, "r" = r, "kmeans.centers" = kmeans.centers, "kmeans.clusters" = kmeans.clusters)
  return(retList)

}

####--------------------------------------------------------------
## sample.bmm.component.mean:  Sample (scalar) means in a particular dimension
##                             from the posterior of a particular component 
##                             in the Beta mixture model
##
## Inputs:  
##
## mu:  a (scalar) shape parameter for the gamma prior distribution over
##      a particular component's u parameter.
## alpha:  a (scalar) rate (i.e., inverse scale) parameter for the gamma prior 
##         distribution over a particular component's u parameter.
## nu:  a (scalar) shape parameter for the gamma prior distribution over
##      a particular component's v parameter.
## beta:  a (scalar) rate (i.e., inverse scale) parameter for the gamma prior 
##         distribution over a particular component's v parameter.
## num.samples:  the number of samples to return
##
## Outputs:
## 
## a vector of length num.samples holding means sampled from the posterior
## distribution of the component and dimension of the Beta mixture model
## indicated by the input parameters.

sample.bmm.component.mean <- function(mu, alpha, nu, beta, num.samples) 
{

  # NB:  In Ma's notation, a gamma is parameterized by a shape
  # parameter mu (or nu) and a rate or inverse scale parameter alpha (or beta).
  # Make sure to use the appropriate parameterization in R.

  # Sample u and v (the Beta parameters) from the Gamma priors
  samples.u <- rgamma(num.samples, shape=mu, rate=alpha)
  samples.v <- rgamma(num.samples, shape=nu, rate=beta)

  # Now the expected means are just u / ( u + v )
  samples.means <- samples.u / ( samples.u + samples.v )

  return(samples.means)

} # End sample.bmm.component.mean

####--------------------------------------------------------------
## sample.bmm.component.proportion:  Sample (scalar) proportions in a 
##                                   particular dimension from the 
##                                   posterior of a particular component in 
##                                   the Beta mixture model
##
## Inputs:  
##
## mu:  a (scalar) shape parameter for the gamma prior distribution over
##      a particular component's u parameter.
## alpha:  a (scalar) rate (i.e., inverse scale) parameter for the gamma prior 
##         distribution over a particular component's u parameter.
## nu:  a (scalar) shape parameter for the gamma prior distribution over
##      a particular component's v parameter.
## beta:  a (scalar) rate (i.e., inverse scale) parameter for the gamma prior 
##         distribution over a particular component's v parameter.
## num.samples:  the number of samples to return
##
## Outputs:
## 
## a vector of length num.samples holding proportions sampled from the 
## posterior distribution of the component and dimension of the Beta 
## mixture model indicated by the input parameters.

sample.bmm.component.proportion <- function(mu, alpha, nu, beta, num.samples) 
{

  # NB:  In Ma's notation, a gamma is parameterized by a shape
  # parameter mu (or nu) and a rate or inverse scale parameter alpha (or beta).
  # Make sure to use the appropriate parameterization in R.

  # Sample u and v (the Beta parameters) from the Gamma priors
  samples.u <- rgamma(num.samples, shape=mu, rate=alpha)
  samples.v <- rgamma(num.samples, shape=nu, rate=beta)

  # Now sample a Beta distribution using these sampled parameters
  samples.proportions <- apply(cbind(samples.u, samples.v), 1, function(row) rbeta(1, shape1=row[1], shape2=row[2]))

  return(samples.proportions)

} # End sample.bmm.component.proportion

####--------------------------------------------------------------
## define.narrowest.interval.including.pt:  Determine the narrowest
##                                          interval within a set of scalars
##                                          that includes a particular scalar
##                                          and a fixed percent of the set.
##
## Inputs:
##
## xs:  a vector of scalars
## pt:  a vector that must be included in the interval
## percentage:  the percent of the scalars that should be included in the
##              interval
##
## Outputs:
##
## a 2-tuple whose first component is the lower bound of the interval
## and whose second component is its upper bound.

define.narrowest.interval.including.pt <- function(xs, pt, percentage){
  width.in.pts <- floor(percentage*length(xs))
  narrowest.width <- Inf
  narrowest.lb <- 0
  narrowest.ub <- 0
  for(i in 1:(length(xs)-width.in.pts)){
    lb <- xs[i]
    if((i+width.in.pts-1)<length(xs)){
      ub <- xs[i+width.in.pts-1]
      if((lb > pt) | (ub<pt)) { next }
      if((ub-lb) < narrowest.width){
        narrowest.width <- (ub-lb)
        narrowest.lb <- lb
        narrowest.ub <- ub
      }
    }
  }
  return(c(narrowest.lb,narrowest.ub))
}

####--------------------------------------------------------------
## bmm.narrowest.mean.interval.about.centers:  Return the narrowest (Bayesian
##    credible) interval of the means that include the centers and
##    a fixed proportion of the probability distribution of the means
##    (e.g., .68 or .95).
##
## Inputs:  
##
## mu:  a D x N.c matrix holding the shape parameters for the gamma prior 
##      distributions over the u parameters.  i.e., mu[d,n] is the 
##      shape parameter governing u[d,n].
##      Introduced in eqn (15).
## alpha:  a D x N.c matrix holding the rate (i.e., inverse scale) 
##         parameters for the gamma prior distributions over the u 
##         parameters.  i.e., mu[d,n] is the rate parameter governing u[d,n].
##         Introduced in eqn (15).
## nu:  a D x N.c matrix holding the hyperparameter values of the 
##      shape parameters for the gamma prior distributions over the 
##      v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
##      Introduced in eqn (16).
## beta:  a D x N.c matrix holding the hyperparameter values of the 
##        rate (i.e., inverse scale) parameters for the gamma prior 
##        distributions over the v parameters.  i.e., beta[d,n] is the rate
##        parameter governing v[d,n]. Introduced in eqn (16).
## proportion:  the percentage of the distribution of each mean to include
##              in the interval (e.g., .68 or .95)
##
## Outputs:
## 
## a list with values lb and ub, where each are D x N.c
## matrices holding the respective lower or upper bounds of the intervals,
## and centers, which is a D x N.c matrix holding the empirical means sampled
## from the posterior distribution

bmm.narrowest.mean.interval.about.centers <- function(mu, alpha, nu, beta, proportion) 
{

  D <- dim(mu)[1]
  N.c <- dim(mu)[2]
      
  num.samples <- 1000

  lb <- matrix(data = 0, nrow=D, ncol=N.c)
  ub <- matrix(data = 0, nrow=D, ncol=N.c)
  centers <- matrix(data = 0, nrow=D, ncol=N.c)

  for(i in 1:N.c){
    for(l in 1:D){
      samples.proportions <- sample.bmm.component.mean(mu[l,i], alpha[l,i], nu[l,i], beta[l,i], num.samples)

      samples.proportions <- sort(samples.proportions, decreasing=FALSE)
  
      centers[l,i] <- mean(samples.proportions)

      bounds <- define.narrowest.interval.including.pt(samples.proportions, centers[l,i], proportion)
      lb[l,i] <- bounds[1]
      ub[l,i] <- bounds[2]
    }
  }
  
  retList <- list("lb" = lb, "ub" = ub, "centers" = centers)
  return(retList)

} # End bmm.narrowest.mean.interval.about.centers 

####--------------------------------------------------------------
## bmm.narrowest.proportion.interval.about.centers:  Return the narrowest 
##    interval of proportions from the posterior distribution that
##    includes their centers a fixed proportion of the posterior distribution
##    (e.g., .68 or .95).
##
## Inputs:  
##
## mu:  a D x N.c matrix holding the shape parameters for the gamma prior 
##      distributions over the u parameters.  i.e., mu[d,n] is the 
##      shape parameter governing u[d,n].
##      Introduced in eqn (15).
## alpha:  a D x N.c matrix holding the rate (i.e., inverse scale) 
##         parameters for the gamma prior distributions over the u 
##         parameters.  i.e., mu[d,n] is the rate parameter governing u[d,n].
##         Introduced in eqn (15).
## nu:  a D x N.c matrix holding the hyperparameter values of the 
##      shape parameters for the gamma prior distributions over the 
##      v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
##      Introduced in eqn (16).
## beta:  a D x N.c matrix holding the hyperparameter values of the 
##        rate (i.e., inverse scale) parameters for the gamma prior 
##        distributions over the v parameters.  i.e., beta[d,n] is the rate
##        parameter governing v[d,n]. Introduced in eqn (16).
## proportion:  the percentage of the distribution of each mean to include
##              in the interval (e.g., .68 or .95)
##
## Outputs:
## 
## a list with values lb and ub, where each are D x N.c
## matrices holding the respective lower or upper bounds of the intervals,
## and centers, which is a D x N.c matrix holding the empirical means sampled
## from the posterior distribution

bmm.narrowest.proportion.interval.about.centers <- function(mu, alpha, nu, beta, proportion) 
{

  D <- dim(mu)[1]
  N.c <- dim(mu)[2]
      
  num.samples <- 1000

  lb <- matrix(data = 0, nrow=D, ncol=N.c)
  ub <- matrix(data = 0, nrow=D, ncol=N.c)
  centers <- matrix(data = 0, nrow=D, ncol=N.c)

  for(i in 1:N.c){
    for(l in 1:D){
      samples.proportions <- sample.bmm.component.proportion(mu[l,i], alpha[l,i], nu[l,i], beta[l,i], num.samples)

      samples.proportions <- sort(samples.proportions, decreasing=FALSE)
  
      centers[l,i] <- mean(samples.proportions)

      bounds <- define.narrowest.interval.including.pt(samples.proportions, centers[l,i], proportion)
      lb[l,i] <- bounds[1]
      ub[l,i] <- bounds[2]
    }
  }
  
  retList <- list("lb" = lb, "ub" = ub, "centers" = centers)
  return(retList)

} # End bmm.narrowest.proportion.interval.about.centers 

####--------------------------------------------------------------
## bmm.component.posterior.predictive.density: Calculate the posterior
##   predictive density in a single dimension for a single component.
##
## Inputs:
##
## x:  the scalar at which to evaluate the posterior predictive density
## mu:  a (scalar) shape parameter for the gamma prior distribution over
##      a particular component's u parameter.
## alpha:  a (scalar) rate (i.e., inverse scale) parameter for the gamma prior 
##         distribution over a particular component's u parameter.
## nu:  a (scalar) shape parameter for the gamma prior distribution over
##      a particular component's v parameter.
## beta:  a (scalar) rate (i.e., inverse scale) parameter for the gamma prior 
##         distribution over a particular component's v parameter.
## pi:  a (scalar) mixing coefficient giving the weight of the component
## num.samples:  the number of samples to use in performing the numerical
##               evaluation of the predictive density
##
## Outputs:
##
## the posterior preditive density at x for the specified component

bmm.component.posterior.predictive.density <- function(x, mu, alpha, nu, beta, pi, num.samples = 1000)
{
  # Posterior predictive density is product of a beta and two gammas
  # Do this numerically.
  
  # Sample parameters from (posterior) gammas
  samples.u <- rgamma(num.samples, shape=mu, rate=alpha)
  samples.v <- rgamma(num.samples, shape=nu, rate=beta)
  
  # Evaluate posterior probability at x for each of these sampled
  # parameters.
  y <- pi * sum(apply(cbind(samples.u, samples.v), 1, function(row) dbeta(x, shape1=row[1], shape2=row[2])))
  
  return(y)
}

####--------------------------------------------------------------
## bmm.posterior.predictive.density: Calculate the posterior
##   predictive density in a single dimension across all components
##
## Inputs:
##
## x:  the scalar at which to evaluate the posterior predictive density
## mu:  an N.c vector holding the shape parameter for each of the N.c 
##      components, which governs the gamma prior distribution over that 
##      component's u parameter.
## alpha:  an N.c vector holding the rate (i.e., inverse scale) parameter 
##         for each of the N.c components, which governs the gamma prior 
##         distribution over that component's u parameter.
## nu:  an N.c vector holding the shape parameter for each of the N.c 
##      components, which governs the gamma prior distribution over that 
##      component's v parameter.
## beta:  an N.c vector holding the rate (i.e., inverse scale) parameter 
##         for each of the N.c components, which governs the gamma prior 
##         distribution over that component's v parameter.
## pi:  an N.c vector holding the mixing coefficients giving the weight 
##      of each of the N.c components
## num.samples:  the number of samples to use in performing the numerical
##               evaluation of the predictive density
## Outputs:
##
## the posterior preditive density at x summed across all components

bmm.posterior.predictive.density <- function(x, mu, alpha, nu, beta, pi, num.samples = 1000)
{
  N.c <- length(mu)
  y <- 0
  for(k in 1:N.c) {
    y <- y + bmm.component.posterior.predictive.density(x, mu[k], alpha[k], nu[k], beta[k], pi[k], num.samples)
  }
  return(y)
}

####--------------------------------------------------------------
## generate.proportion.data.set:  Generate a data set with clusters of 
##    Beta-derived proportions
##
## Generate N.c clusters of proportions, where each dimension is
## independently generated by a beta distribution.
##
## Inputs:
##
## N.c:  the number of clusters
## N:  the total number of points (which will roughly proportionally
##     allocated across the N.c clusters)
## D:  the number of dimensions of the data set.
##
## Outputs:
##
## an N x D matrix of proportions
##

generate.proportion.data.set <- function(N.c, N, D)
{
  m <- matrix(data = 0, nrow=N, ncol=D)
  num.items.per.cluster <- floor(N/N.c)
  nxt.indx <- 1
  u <- 0
  v <- 0
  for(i in 1:N.c) {
    for(j in 1:D) {
      # Generate a random center for this cluster in this dimension
      u <- runif(1, min=10, max=1000)
      v <- i*runif(1, min=10, max=1000)
      m[nxt.indx:(nxt.indx+num.items.per.cluster-1),j] <- rbeta(num.items.per.cluster, u, v)
    }
    nxt.indx <- nxt.indx + num.items.per.cluster
  }

  # If N does not evenly divide N.c, generate the remaining points
  if(nxt.indx <= N) {
    for(j in 1:D) {
      m[nxt.indx:(nxt.indx+(N-nxt.indx)),j] <- rbeta(N+1-nxt.indx, u, v)
    }
  }

  # Keep proportions away from 0 or 1
  delta <- .Machine$double.eps
  m[m > ( 1 - delta )] <- 1 - delta
  m[m < delta] <- delta

  return(m)
}


####--------------------------------------------------------------
## bmm.plot.1d:  Overlay the posterior predictive density on a histogram
##               of the data.
##
## Inputs:
## X:  an N x D matrix with rows being the items to cluster.
##     All entries are assumed to be proportions (i.e., between 0 and 1).
##     Notice that there are no summation restrictions--i.e., proportions
##     do not sum to unity across an item's dimensions.
## mu:  a D x N.c matrix holding the values of the 
##      shape parameters for the gamma prior distributions over the 
##      u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
##      Introduced in eqn (15).
## alpha:  a D x N.c matrix holding the values of the 
##         rate (i.e., inverse scale) parameters for the gamma prior 
##         distributions over the u parameters.  i.e., mu[d,n] is the rate
##         parameter governing u[d,n]. Introduced in eqn (15).
## nu:  a D x N.c matrix holding the values of the 
##      shape parameters for the gamma prior distributions over the 
##      v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
##      Introduced in eqn (16).
## beta:  a D x N.c matrix holding the values of the 
##        rate (i.e., inverse scale) parameters for the gamma prior 
##        distributions over the v parameters.  i.e., beta[d,n] is the rate
##        parameter governing v[d,n]. Introduced in eqn (16).
## E.pi:  the D-vector holding the values E[pi], i.e., the expected values
##        of the mixing coefficients, defined in eqn (53).
## r:  the N x N.c matrix of responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## title:  plot title
## xlab:  x label
## ylab:  y label
##
## Outputs:
##
## ggplot object overlaying the posterior predictive density on a histogram
## of the data.

bmm.plot.1d <- function(X, mu, alpha, nu, beta, E.pi, r, title, xlab, ylab)
{
  bin.width <- .025
  N.c <- dim(mu)[2]

  proportions <- data.frame(x=X[,1], row.names=NULL, stringsAsFactors=NULL)
  
  # Generate (x,y) values of the posterior predictive density
  n <- 1000
  y <- rep.int(0, n)
  # Don't evaluate at x=0 or x=1, which will blow up
  x <- seq(1/n, 1-(1/n), length=n)
  ym <- matrix(data=0,nrow=N.c,ncol=n)
  num.iterations <- 1000
  for (k in 1:N.c) {
    for (i in 1:n) {

      # Evaluate posterior probability at x.
      ym[k,i] <- bmm.component.posterior.predictive.density(x[i], mu[1,k], alpha[1,k], nu[1,k], beta[1,k], E.pi[k], num.samples = num.iterations)
  
      y[i] <- y[i] + ym[k,i]
    }
  }
  
  max.posterior.density <- max(ym)
  # Overlay the histogram of the data on to the beta mixture model fit.

  # Set max posterior to max of splinefun.
  limits <- data.frame(x=c(min(x),max(x)))

  f <- splinefun(x,y)
  max.posterior.density <- max(unlist(lapply(seq(from=limits$x[1],to=limits$x[2],by=10^-3), f)))

  g <- ggplot(data = proportions, aes(x)) + ggtitle(title) + xlab(xlab) + ylab(ylab)

  g <- g + theme_bw() + theme(axis.line = theme_segment(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) 

  num.breaks <- ceiling(1/bin.width) + 1
  breaks <- unlist(lapply(0:(num.breaks-1), function(x) x/(num.breaks-1)))
  g <- g + geom_histogram(data=proportions, mapping=aes(x), fill="white", colour="black", breaks=breaks)

  tmp <- print(g)
  max.density <- max(tmp[["data"]][[1]]$ymax)

  scale <- (max.density/max.posterior.density)

  vertical.offset <- .1 * max.posterior.density

  for(k in 1:N.c) {
      
    f <- splinefun(x,scale*ym[k,])

    g <- g + stat_function(data = limits, fun=f, mapping=aes(x), lty="dashed")

  }

  f <- splinefun(x,scale*y)

  g <- g + stat_function(data = limits, fun=f, mapping=aes(x))
  
  xmin <- -0.05
  xmax <-  1.05

  g <- g + coord_cartesian(ylim=c(0, max.density*1.1), xlim=c(xmin,xmax))

  return(g)
} # End bmm.plot.1d

my.ellipse <- function(hlaxa = 1, hlaxb = 1, theta = 0, xc = 0, yc = 0, newplot = F,npoints = 100, ...){
  a <- seq(0, 2 * pi, length = npoints + 1)
  x <- hlaxa * cos(a)
  y <- hlaxb * sin(a)
  alpha <- angle(x, y)
  rad <- sqrt(x^2 + y^2)
  xp <- rad * cos(alpha + theta) + xc
  yp <- rad * sin(alpha + theta) + yc
  return (cbind(data.frame(x=xp, y=yp)))
}
  
angle <- function (x, y){
  angle2 <- function(xy) {
    x <- xy[1]
    y <- xy[2]
    if (x > 0) {
      atan(y/x)
    }else {
      if (x < 0 & y != 0) {
        atan(y/x) + sign(y) * pi
      }else {
        if (x < 0 & y == 0) {
          pi
        }else {
          if (y != 0) {
            (sign(y) * pi)/2
          }else {
            NA
          }
        }
      }
    }
  }
  apply(cbind(x, y), 1, angle2)
}

####--------------------------------------------------------------
## bmm.plot.2d:  Plot data highlighted according to clustering
##
## Plot the data with each item's respective color/symbol indicating
## its cluster.  Also display "1-sigma" contour intervals for the
## standard error of the mean for the standard error.
##
## Inputs:
## X:  an N x D matrix with rows being the items to cluster.
##     All entries are assumed to be proportions (i.e., between 0 and 1).
##     Notice that there are no summation restrictions--i.e., proportions
##     do not sum to unity across an item's dimensions.
## mu:  a D x N.c matrix holding the values of the 
##      shape parameters for the gamma prior distributions over the 
##      u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
##      Introduced in eqn (15).
## alpha:  a D x N.c matrix holding the values of the 
##         rate (i.e., inverse scale) parameters for the gamma prior 
##         distributions over the u parameters.  i.e., mu[d,n] is the rate
##         parameter governing u[d,n]. Introduced in eqn (15).
## nu:  a D x N.c matrix holding the values of the 
##      shape parameters for the gamma prior distributions over the 
##      v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
##      Introduced in eqn (16).
## beta:  a D x N.c matrix holding the values of the 
##        rate (i.e., inverse scale) parameters for the gamma prior 
##        distributions over the v parameters.  i.e., beta[d,n] is the rate
##        parameter governing v[d,n]. Introduced in eqn (16).
## E.pi:  the D-vector holding the values E[pi], i.e., the expected values
##        of the mixing coefficients, defined in eqn (53).
## r:  the N x N.c matrix of responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## title:  plot title
## xlab:  x label
## ylab:  y label
##
## Outputs:
##
## ggplot object displaying the clustered data and "1-sigma" contour intervals
## for the standard error of the mean and for the standard error

bmm.plot.2d <- function(X, mu, alpha, nu, beta, E.pi, r, title, xlab, ylab)
{
  N <- dim(X)[1]
  N.c <- dim(mu)[2]

  # width = 1 std dev
  suppressPackageStartupMessages(library("NORMT3")) # for erf
  width <- as.double(erf(1/sqrt(2)))  
 
  width <- sqrt(width)

  # Calculate standard error of the means
  SEM.res <- bmm.narrowest.mean.interval.about.centers(mu, alpha, nu, beta, width)
  SEM.centers <- t(SEM.res$centers)
  SEMs.lb <- t(SEM.res$lb)
  SEMs.ub <- t(SEM.res$ub)

  # Calculate standard errors
  std.dev.res <- bmm.narrowest.proportion.interval.about.centers(mu, alpha, nu, beta, width)
  std.dev.centers <- t(std.dev.res$centers)
  std.dev.lb <- t(std.dev.res$lb)
  std.dev.ub <- t(std.dev.res$ub)

  # Make sure none of the clusters = 0 or 1, which would be used for black.
  # Let's reserve that for highlighting.

  clusters <- rep(0, N)
  for(n in 1:N) {
    max.cluster <- 0
    max.assignment <- -1
    for(k in 1:N.c) {
      if ( r[n,k] > max.assignment ) {
        max.assignment <- r[n,k]
        # + 2 ensures none of the clusters has number 0 or 1,
        # which would make it black when cluster is used as colour
        max.cluster <- ( k + 2 )
      }
    }
    clusters[n] <- max.cluster
  }

  proportions <- data.frame(x=X[,1], y=X[,2], row.names=NULL, stringsAsFactors=NULL)

  centers <- data.frame(x=SEM.centers[,1], y=SEM.centers[,2], row.names=NULL, stringsAsFactors=NULL)

  # R CMD check does not understand that the right-hand side of the = in
  # aes commands is referring, not to global or local variables, but to
  # verbatim names in the data argument.  Get around this by creating
  # variables with those names--which is irrelevant to the ggplot code.
  x <- NULL
  y <- NULL
  g <- ggplot(data = proportions, aes(x=x, y=y)) + ggtitle(title) + xlab(xlab) + ylab(ylab) + geom_point(data = proportions, aes(x=x, y=y), shape=clusters, colour=clusters) + geom_point(data = centers, aes(x=x, y=y), size=3, colour='blue') 

  g <- g + theme_bw() + theme(axis.line = theme_segment(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) 
  
  # Add the "1-sigma" contour intervals for the standard error _of the mean_
  for(k in 1:N.c) {
    # i.e., provide the 1-sigma confidence interval of q(mu_k) =
    # \int q(mu_k, lambda_k) dlambda_k.
   
    xc <- SEMs.lb[k,1] + ((SEMs.ub[k,1] - SEMs.lb[k,1])/2)
    yc <- SEMs.lb[k,2] + ((SEMs.ub[k,2] - SEMs.lb[k,2])/2)
  
    ell <- my.ellipse(hlaxa = ((SEMs.ub[k,1] - SEMs.lb[k,1])/2), hlaxb = ((SEMs.ub[k,2] - SEMs.lb[k,2])/2), xc = xc, yc = yc)
  
    df_ell <- cbind(data.frame(x=ell$x, y=ell$y))
  
    g <- g + geom_path(data=df_ell, aes(x=x,y=y), colour="blue", linetype=2)
  }
  
  # Add the "1-sigma" contour intervals for the standard error  
  for(k in 1:N.c) {
    # i.e., provide the 1-sigma confidence interval of the posterior
    # predictive density (Bishop eqn 10.81)--a Student t.
  
    xc <- std.dev.lb[k,1] + ((std.dev.ub[k,1] - std.dev.lb[k,1])/2)
    yc <- std.dev.lb[k,2] + ((std.dev.ub[k,2] - std.dev.lb[k,2])/2)
  
    ell <- my.ellipse(hlaxa = ((std.dev.ub[k,1] - std.dev.lb[k,1])/2), hlaxb = ((std.dev.ub[k,2] - std.dev.lb[k,2])/2), xc = xc, yc = yc)
  
    df_ell <- cbind(data.frame(x=ell$x, y=ell$y))
  
    g <- g + geom_path(data=df_ell, aes(x=x,y=y))
  }  

  xmin <- -0.05
  xmax <-  1.05

  g <- g + coord_cartesian(xlim=c(xmin,xmax), ylim=c(xmin,xmax))

  return(g)

} # End bmm.plot.2d

## Begin binomial mixture model

####--------------------------------------------------------------
## binomial.bmm.fixed.num.components: Use a variational Bayesian approach to fit 
## a mixture of binomial distributions to proportion data, without dropping
## any components/clusters.  To instead automatically determine the number of
## components, use binomial.bmm, which invokes this function.
##
## Inputs:
## X:  an N x D matrix with rows holding the number of successes
##              in each dimension.
##              All entries are assumed to be counts--not proportions.
## eta:  an N x D matrix with rows holding the number of trials
##                in each dimension.
##                All entries are assumed to be counts--not proportions.
## N.c:  the number of components/clusters to attempt
## r:  the N x N.c matrix of initial responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## a:  a N.c x D matrix holding the _initial_ values of the 
##     shape parameters a for the beta prior distributions over the 
##     mu parameters, i.e., Beta(mu[k,m]; a0[k,m], b0[k,m])  
## b:  a N.c x D matrix holding the _initial_ values of the 
##     shape parameters b for the beta prior distributions over the 
##     mu parameters, i.e., Beta(mu[k,m]; a0[k,m], b0[k,m])  
## alpha:  a vector with D components holding the _initial_ values of the
##         parameters of the Dirichlet distribution over the mixing 
##         coefficients pi.  
## a0, b0, alpha0:  the hyperparameters corresponding to the
##                               above initial values (and with the same
##                               respective matrix/vector dimensionality).
## convergence.threshold:  minimum absolute difference between mixing
##                         coefficient (expected) values across consecutive
##                         iterations to reach converge.
## max.iterations:  maximum number of iterations to attempt
## verbose:  output progress in terms of mixing coefficient (expected) values
##           if 1.
## Outputs: a list with the following entries
## retVal:  0 indicates successful convergence; -1 indicates a failure
##          to converge.
## a:  a N.c x D matrix holding the _converged final_ values of the 
##     shape parameters a for the beta posterior distributions over the 
##     mu parameters, i.e., Beta(mu[k,m]; a[k,m], b[k,m])  
## b:  a N.c x D matrix holding the _converged final_ values of the 
##     shape parameters b for the beta posterior distributions over the 
##     mu parameters, i.e., Beta(mu[k,m]; a[k,m], b[k,m])  
## alpha:  a vector with D components holding the _converged final_ values of 
##         the parameters of the posterior Dirichlet distribution over the 
##         mixing coefficients pi.  
## r:  the N x N.c matrix of responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## num.iterations:  the number of iterations required to reach convergence.
## ln.rho:  an N x N.c matrix holding the ln[rho], as defined in eqn (32).
## E.lnpi:  the D-vector holding the values E[ln pi], defined following 
##          eqn (51).
## E.pi:  the D-vector holding the values E[pi], i.e., the expected values
##        of the mixing coefficients, defined in eqn (53).

binomial.bmm.fixed.num.components <- function(X, eta, N.c, r, a, b, alpha, a0, b0, alpha0, convergence.threshold = 10^-4, max.iterations = 10000, verbose = 0)
{
  N <- dim(X)[1]
  D <- dim(X)[2]

  E.pi.prev <- rep(0, N.c)
  
  iteration <- 0

  # Apply variational bayesian approach to binomial mixture modeling
  # until convergence

  Nk <- colSums(r)
  if (any(is.na(Nk))) {
    print(Nk)
    stop("Nk is NA")
  }
    
  # Calculate xbar_km
  Nk.xbar <- t(r) %*% X
  if (any(is.na(Nk.xbar))) {
    print(Nk.xbar)
    print("r")
    print(r)
    print("X")
    print(X)
    stop("Nk.xbar is NA")
  }
    
  # Calculate etabar_km
  Nk.etabar <- t(r) %*% eta
  if (any(is.na(Nk.etabar))) {
    print(Nk.etabar)
    stop("Nk.etabar is NA")
  }

  lb.prev <- -Inf

  while(TRUE) {  

    # Loop:
    #    (1) Update vars (a, b, alpha)

    # Calculate alpha_k
    alpha <- alpha0 + Nk
    if (any(is.na(alpha))) {
      print(alpha)
      stop("alpha is NA")
    }
    
    # Calculate a_km
    # Calculate b_km
    a <- ( a0 - 1 ) + Nk.xbar
    b <- ( b0 - 1 ) + Nk.etabar - Nk.xbar
    if (any(is.na(a))) {
      print(a)
      stop("a is NA")
    }
   
    if (any(is.na(b))) {
      print(b)
      stop("b is NA")
    } 

    if ( any(a <= 0) ) {
      print(a)
      #print(Nk.xbar)
      #print(r)
      print(alpha/sum(alpha))
      stop("a <= 0")
    }

    if ( any(b <= 0) ) {
      stop("b <= 0")
    }

    if ( any(alpha < 0) ) {
      stop("alpha < 0")
    }    
    
    #    (2) Compute expectations
  
    # Calculate E_mu[ln mu_km]
    E.ln.mu <- digamma(a) - digamma(a + b)
    if (any(is.na(E.ln.mu))) {
      print(E.ln.mu)
      print("a")
      print(a)
      print("b")
      print(b)
      stop("E.ln.mu is NA")
    }
  
    # Calculate E_mu[ln ( 1 - mu_km )]
    E.ln.one.minus.mu <- digamma(b) - digamma(a + b)
    if (any(is.na(E.ln.one.minus.mu))) {
      print(E.ln.one.minus.mu)
      stop("E.ln.one.minus.mu is NA")
    }
    
    # Calculate E_pi[ln pi_k]
    E.lnpi <- digamma(alpha) - digamma(sum(alpha))
    if (any(is.na(E.lnpi))) {
      print(E.lnpi)
      stop("E.lnpi is NA")
    }

    
    #    (3) Compute responsibilities
    
    # Calculate responsibilities, E_Z[z_nk]


    one.matrix <- matrix(data=1, nrow=D, ncol=N.c)
    
    ln.rho <- matrix(data=E.lnpi, nrow=N, ncol=N.c, byrow=TRUE)
    ln.rho <- ln.rho + X %*% t(E.ln.mu) + eta %*% t(E.ln.one.minus.mu) - X %*% t(E.ln.one.minus.mu) 
    ln.rho <- ln.rho + ( lgamma(eta + 1) %*% one.matrix ) - ( lgamma(eta - X + 1) %*% one.matrix ) - ( lgamma(X + 1) %*% one.matrix )
    if (any(is.na(ln.rho))) {
      print(ln.rho)
      stop("ln.rho is NA")
    }

    r <- matrix(data = 0, nrow=N, ncol=N.c)
    for(n in 1:N) {
      if(any(is.na(ln.rho[n,]))) {
        r[n,] <- rep(NA, N.c)
        next
      }
      row.sum <- log(sum(exp(ln.rho[n,] - max(ln.rho[n,])))) + max(ln.rho[n,], na.rm=TRUE)
      for(k in 1:N.c) { r[n,k] = exp(ln.rho[n,k] - row.sum) }
    }

    if (any(is.na(r))) {
      print(r)
      stop("r is NA")
    }
    
    #    (4) Compute statistics

  
    # Calculate N_k
    r <- r + 10^-9

    Nk <- colSums(r)
    if (any(is.na(Nk))) {
      print(Nk)
      stop("Nk is NA")
    }
    
    # Calculate xbar_km
    Nk.xbar <- t(r) %*% X
    if (any(is.na(Nk.xbar))) {
      print(Nk.xbar)
      print("r")
      print(r)
      print("X")
      print(X)
      stop("Nk.xbar is NA")
    }
    
    # Calculate etabar_km
    Nk.etabar <- t(r) %*% eta
    if (any(is.na(Nk.etabar))) {
      print(Nk.etabar)
      stop("Nk.etabar is NA")
    }
    
    #    (5) Compute bound
    
    # Compute variational lower bound
    tmp.mat <- lgamma(eta + 1) - lgamma(X + 1) - lgamma(eta - X + 1)
    # NB:  last two terms are element-wise multiplication by r, not
    # matrix-matrix multiplication
    E.ln.p.X.Z.mu <- sum(t(tmp.mat) %*% r) + sum((X %*% t(E.ln.mu)) * r) + sum(((eta - X) %*% t(E.ln.one.minus.mu)) * r)

    E.ln.p.Z.pi <- sum(r %*% E.lnpi)

    E.ln.p.pi <- lgamma(sum(alpha0)) - sum(lgamma(alpha0)) + sum( (alpha0 - 1) * E.lnpi )
    
    E.ln.p.mu <- sum( (lgamma(a0 + b0)) - (lgamma(a0)) - (lgamma(b0)) + ((a0 - 1) * E.ln.mu) + ((b0-1) * E.ln.one.minus.mu) )

    E.ln.q.Z <- 0
    for(n in 1:N) {
      for(k in 1:N.c) {
        if ( r[n,k] != 0 ) {
          E.ln.q.Z <- E.ln.q.Z + r[n,k] * log(r[n,k])
        }
      }
    }

    E.ln.q.pi <- lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1)*E.lnpi)
    
    E.ln.q.mu <- sum( lgamma(a+b) - lgamma(a) - lgamma(b) + ((a-1)*E.ln.mu) + ((b-1)*E.ln.one.minus.mu) )

    if (any(is.na(E.ln.p.X.Z.mu))) {
      stop("E.ln.p.X.Z.mu is NA")
    }

    if (any(is.na(E.ln.p.Z.pi))) {
      stop("E.ln.p.Z.pi is NA")
    }

    if (any(is.na(E.ln.p.pi))) {
      stop("E.ln.p.pi is NA")
    }

    if (any(is.na(E.ln.p.mu))) {
      stop("E.ln.p.mu is NA")
    }

    if (any(is.na(E.ln.q.Z))) {
      stop("E.ln.q.Z is NA")
    }

    if (any(is.na(E.ln.q.pi))) {
      stop("E.ln.q.pi is NA")
    }

    if (any(is.na(E.ln.q.mu))) {
      stop("E.ln.q.mu is NA")
    }

    lb <- E.ln.p.X.Z.mu + E.ln.p.Z.pi + E.ln.p.pi + E.ln.p.mu - E.ln.q.Z - E.ln.q.pi - E.ln.q.mu 
    
    iteration <- iteration + 1
    # if ( iteration > 100 ) { break }

    # cat(lb, "\n")

    E.pi <- alpha / sum(alpha)
    if(verbose) {
      cat("lb = ", lb, " pi = ", E.pi, "\n")
    }

    if ( lb.prev > lb ) {
      cat(sprintf("lb decreased from %f to %f!\n", lb.prev, lb))
    }
    if ( abs( ( lb - lb.prev ) / lb ) < convergence.threshold ) { break }

    lb.prev <- lb

  } # End inner while(TRUE)

  retList <- list("retVal" = 0, "a" = a, "b" = b, "alpha" = alpha, "r" = r, "num.iterations" = iteration, "ln.rho" = ln.rho, "E.lnpi" = E.lnpi, "E.pi" = E.pi)

  return(retList)
} # End binomial.bmm.fixed.num.components function

####--------------------------------------------------------------
## binomial.bmm: Use a variational Bayesian approach to fit a mixture of 
## binomial distributions to count data.
##
## NB: Clustering is first run to convergence using 
## binomial.bmm.fixed.num.components.
## If any components have probability less than pi.threshold, they are
## discarded and the clustering is run to convergence again.
##
##
## Inputs:
## successes:  an N x D matrix with rows holding the number of successes
##              in each dimension.
##              All entries are assumed to be counts--not proportions.
## total.trials:  an N x D matrix with rows holding the number of trials
##                in each dimension.
##                All entries are assumed to be counts--not proportions.
## N.c:  the number of components/clusters to attempt
## r:  the N x N.c matrix of initial responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## mu:  a D x N.c matrix holding the _initial_ values of the 
##      shape parameters for the gamma prior distributions over the 
##      u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
##      NB:  this is the initial value mu, which is updated upon iteration.
##      It is not (necessarily) the same as the hyperparameter mu0, which
##      is unchanged by iteration.
##      Introduced in eqn (15).
## alpha:  a D x N.c matrix holding the _initial_ values of the 
##         rate (i.e., inverse scale) parameters for the gamma prior 
##         distributions over the u parameters.  i.e., mu[d,n] is the rate
##         parameter governing u[d,n]. Introduced in eqn (15).
##         NB:  this is the initial value alpha, which is updated upon iteration.
##         It is not (necessarily) the same as the hyperparameter alpha0, which
##         is unchanged by iteration.
## nu:  a D x N.c matrix holding the _initial_ values of the 
##      shape parameters for the gamma prior distributions over the 
##      v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
##      Introduced in eqn (16).
##      NB:  this is the initial value nu, which is updated upon iteration.
##      It is not (necessarily) the same as the hyperparameter nu0, which
##      is unchanged by iteration.
## beta:  a D x N.c matrix holding the _initial_ values of the 
##        rate (i.e., inverse scale) parameters for the gamma prior 
##        distributions over the v parameters.  i.e., beta[d,n] is the rate
##        parameter governing v[d,n]. Introduced in eqn (16).
##        NB:  this is the initial value beta, which is updated upon iteration.
##        It is not (necessarily) the same as the hyperparameter beta0, which
##        is unchanged by iteration.
## c:  a vector with D components holding the _initial_ values of the
##     parameters of the Dirichlet distribution over the mixing 
##     coefficients pi.  Introduced in eqn (19).
##     NB: this is the initial value c, which is updated upon iteration.
##     It is not (necessarily) the same as the hyperparameter c0, which 
##     is unchanged by iteration.
## mu0, alpha0, nu0, beta0, c0:  the hyperparameters corresponding to the
##                               above initial values (and with the same
##                               respective matrix/vector dimensionality).
## convergence.threshold:  minimum absolute difference between mixing
##                         coefficient (expected) values across consecutive
##                         iterations to reach converge.
## max.iterations:  maximum number of iterations to attempt
## verbose:  output progress in terms of mixing coefficient (expected) values
##           if 1.
## pi.threshold:  discard any cluster with weight/mixing coefficient less
##                than pi.threshold _following_ convergence.  
## Outputs: a list with the following entries
## retVal:  0 indicates successful convergence; -1 indicates a failure
##          to converge.
## mu:  a D x N.c matrix holding the _converged final_ values of the 
##      shape parameters for the gamma prior distributions over the 
##      u parameters.  i.e., mu[d,n] is the shape parameter governing u[d,n].
##      Introduced in eqn (15).
## alpha:  a D x N.c matrix holding the _converged final_ values of the 
##         rate (i.e., inverse scale) parameters for the gamma prior 
##         distributions over the u parameters.  i.e., mu[d,n] is the rate
##         parameter governing u[d,n]. Introduced in eqn (15).
## nu:  a D x N.c matrix holding the _converged final_ values of the 
##      shape parameters for the gamma prior distributions over the 
##      v parameters.  i.e., nu[d,n] is the shape parameter governing v[d,n].
##      Introduced in eqn (16).
## beta:  a D x N.c matrix holding the _converged final_ values of the 
##        rate (i.e., inverse scale) parameters for the gamma prior 
##        distributions over the v parameters.  i.e., beta[d,n] is the rate
##        parameter governing v[d,n]. Introduced in eqn (16).
## c:  a vector with D components holding the _converged final_ values of the
##     parameters of the Dirichlet distribution over the mixing 
##     coefficients pi.  Introduced in eqn (19).
## r:  the N x N.c matrix of responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## num.iterations:  the number of iterations required to reach convergence.
## ln.rho:  an N x N.c matrix holding the ln[rho], as defined in eqn (32).
## E.lnpi:  the D-vector holding the values E[ln pi], defined following 
##          eqn (51).
## E.pi:  the D-vector holding the values E[pi], i.e., the expected values
##        of the mixing coefficients, defined in eqn (53).

binomial.bmm <- function(successes, total.trials, N.c, r, a, b, alpha, a0, b0, alpha0, convergence.threshold = 10^-4, max.iterations = 10000, verbose = 0, pi.threshold = 10^-2)
{

  total.iterations <- 0 
  D <- dim(successes)[2]

  while(TRUE) {

    if(verbose){
      print(r)
    }
    
    bmm.res <- binomial.bmm.fixed.num.components(successes, total.trials, N.c, r, a, b, alpha, a0, b0, alpha0, convergence.threshold = 10^-4, max.iterations = 10000, verbose = verbose)
    if(bmm.res$retVal != 0) {
      stop("Failed to converge!\n")
    }

    a <- bmm.res$a
    b <- bmm.res$b
    alpha <- bmm.res$alpha
    E.pi <- bmm.res$E.pi

    total.iterations <- total.iterations + bmm.res$num.iterations

    ln.rho <- bmm.res$ln.rho
    E.lnpi <- bmm.res$E.lnpi

    do.inner.iteration <- FALSE

    apply.min.items.condition <- TRUE

    if((apply.min.items.condition == TRUE) & (N.c > 1)) {
      non.zero.indices <- E.pi > pi.threshold
      N = length(successes[,1])
      clusters <- rep(0, N)
      for(n in 1:N) {
        max.cluster <- 0
        max.assignment <- -1
        for(k in 1:N.c) {
          if ( r[n,k] > max.assignment ) {
            max.assignment <- r[n,k]
            max.cluster <- k
          }
        }
        clusters[n] <- max.cluster
      }

      if ( any(non.zero.indices==FALSE) ) {

        do.inner.iteration <- TRUE

        numeric.indices <- (1:N.c)

        E.pi <- E.pi[non.zero.indices]
        E.lnpi <- E.lnpi[non.zero.indices]

        N.c <- length(E.pi)
        r <- matrix(r[,non.zero.indices], nrow=N, ncol=N.c)
        ln.rho <- matrix(ln.rho[,non.zero.indices], nrow=N, ncol=N.c)
        # Need to renormalize r--do it gently.
        for(n in 1:N) {
           if(any(is.na(ln.rho[n,]))) {
             r[n,] <- rep(NA, N.c)
             next
           }
          
           row.sum <- log(sum(exp(ln.rho[n,] - max(ln.rho[n,])))) + max(ln.rho[n,])
           for(k in 1:N.c) { r[n,k] = exp(ln.rho[n,k] - row.sum) }
        }
        a <- matrix(a[non.zero.indices,], nrow=N.c, ncol=D)
        b <- matrix(b[non.zero.indices,], nrow=N.c, ncol=D)
        alpha <- alpha[non.zero.indices,drop=FALSE]
        a0 <- matrix(a0[non.zero.indices,], nrow=N.c, ncol=D)
        b0 <- matrix(b0[non.zero.indices,], nrow=N.c, ncol=D)
        alpha0 <- alpha0[non.zero.indices,drop=FALSE]
      }
    } # End apply.min.items.condition

    if(do.inner.iteration == FALSE) { break }

  }

  retList <- list("retVal" = 0, "a" = a, "b" = b, "alpha" = alpha, "r" = r, "num.iterations" = total.iterations, "ln.rho" = ln.rho, "E.lnpi" = E.lnpi, "E.pi" = E.pi)

  return(retList)

} # End binomial.bmm function


####--------------------------------------------------------------
## init.binomial.bmm.hyperparameters:  Initialize the binomial mixture model
##                                     hyperparameters (to be passed to 
##                                     binomial.bmm)
##
## NB:  This provides what should be a generally reasonable initialization of
##      hyperparameters.  However, better results may be obtained by
##      tuning these in an application-specific manner.
##
## Inputs:
## successess:  an N x D matrix with rows holding the number of successes
##              in each dimension.
##              All entries are assumed to be counts--not proportions.
## total.trials:  an N x D matrix with rows holding the number of trials
##                in each dimension.
##                All entries are assumed to be counts--not proportions.
## N.c:  the number of components/clusters to attempt
##
## Outputs:
## a0:  a N.c x D matrix holding the hyperparameter values of the 
##      shape parameters a for the beta prior distributions over the 
##      mu parameters, i.e., Beta(mu[k,m]; a0[k,m], b0[k,m])  
## b0:  a N.c x D matrix holding the hyperparameter values of the 
##      shape parameters b for the beta prior distributions over the 
##      mu parameters, i.e., Beta(mu[k,m]; a0[k,m], b0[k,m])  
## alpha0:  a vector with D components holding the hyperparameter values of the
##          parameters of the Dirichlet distribution over the mixing 
##          coefficients pi.  
init.binomial.bmm.hyperparameters <- function(successes, total.trials, N.c) 
{
  D <- dim(successes)[2]

  # Choose the initial parameters.  All priors are betas and/or dirichlets.
  # Hence choosing prior parameters = 1 gives a completely flat prior.
  a0 <- matrix(data=1, nrow=N.c, D)
  b0 <- a0

  alpha0 <- rep(0.001, N.c)

  retList <- list("a0" = a0, "b0" = b0, "alpha0" = alpha0)
  return(retList)

} # End init.binomial.bmm.hyperparameters


####--------------------------------------------------------------
## init.binomial.bmm.parameters:  Initialize the binomial mixture model 
##                                parameters (to be passed to binomial.bmm)
##
## Initialize parameters such that expected proportions have the values
## determined by an initial k-means clustering.
##
## NB:  This provides what should be a generally reasonable initialization of
##      hyperparameters.  However, better results may be obtained by
##      tuning these in an application-specific manner.
##
##
## Inputs:
## successess:  an N x D matrix with rows holding the number of successes
##              in each dimension.
##              All entries are assumed to be counts--not proportions.
## total.trials:  an N x D matrix with rows holding the number of trials
##                in each dimension.
##                All entries are assumed to be counts--not proportions.
## N.c:  the number of components/clusters to attempt
## a0:  a N.c x D matrix holding the hyperparameter values of the 
##      shape parameters a for the beta prior distributions over the 
##      mu parameters, i.e., Beta(mu[k,m]; a0[k,m], b0[k,m])  
## b0:  a N.c x D matrix holding the hyperparameter values of the 
##      shape parameters b for the beta prior distributions over the 
##      mu parameters, i.e., Beta(mu[k,m]; a0[k,m], b0[k,m])  
## alpha0:  a vector with D components holding the hyperparameter values of the
##          parameters of the Dirichlet distribution over the mixing 
##          coefficients pi.  
##
## Outputs:
## a:  a N.c x D matrix holding the _initial_ values of the 
##     shape parameters a for the beta prior distributions over the 
##     mu parameters, i.e., Beta(mu[k,m]; a0[k,m], b0[k,m])  
## b:  a N.c x D matrix holding the _initial_ values of the 
##     shape parameters b for the beta prior distributions over the 
##     mu parameters, i.e., Beta(mu[k,m]; a0[k,m], b0[k,m])  
## alpha:  a vector with D components holding the _initial_ values of the
##         parameters of the Dirichlet distribution over the mixing 
##         coefficients pi.  
## r:  the N x N.c matrix of initial responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## kmeans.clusters:  an N-vector giving the assignment of each of the N
##                   items to a cluster, as determined by kmeans.
## kmeans.centers:  an N.c x D matrix holding the centers of the N.c
##                  clusters/components determined by kmeans
init.binomial.bmm.parameters <- function(successes, total.trials, N.c, a0, b0, alpha0, verbose=0) 
{
  N <- dim(successes)[1]
  D <- dim(successes)[2]

  # Initialize the responsibilities r_ni with K-means.
  kmeans.out <- kmeans(successes/total.trials, N.c, nstart=1000)
  kmeans.clusters <- kmeans.out$cluster
  kmeans.centers <- kmeans.out$centers

  # Sort the centers to ease comparison across runs.  No we sort thes
  # clustered results.
  # Besides, be careful!  This upset the relationship between kmeans.clusters
  # and kmeans.centers!  #kmeans.centers <- kmeans.centers[do.call(order, lapply(1:ncol(kmeans.centers), function(i) kmeans.centers[, i])), ]
  if(verbose){
    cat("kmeans initialization:\n")
    write.table(kmeans.centers, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
  }

  r <- matrix(data=0, nrow=N, ncol=N.c)
  for(i in 1:N) {
    r[i,kmeans.clusters[i]]<- 1
  }

  # Set initial a and b (NB: not the prior hyperparameters) according
  # to the means found by kmeans.  Keep mu and nu at their
  # prior values.  To do so:
  # E[proportion] = a / ( a + b )
  # Var[proportion] = a b / [ ( a + b )^2 ( a + b + 1 ) ]

  # Ignore the empirical std dev from kmeans and set the std dev to 0.3
  std.dev <- 0.3
  sigma2 <- std.dev^2

  # These may be solved algebraically to give:
  # a = b * mu / ( 1 - mu ) = b gamma   with gamma = mu / ( 1 - mu )
  # b = ( mu^2 - sigma^2 * gamma ) / [ sigma^2 gamma ( 1 + gamma ) ]
  mu <- kmeans.centers
  gamma <- mu / ( 1 - mu )
  b <- ( mu^2 - sigma2 * gamma ) / ( sigma2 * gamma * (1 + gamma) )
  a <- b * gamma

  alpha <- alpha0 + colSums(r, na.rm=TRUE) 
  
  retList <- list("a" = a, "b" = b, "alpha" = alpha, "r" = r, "kmeans.centers" = kmeans.centers, "kmeans.clusters" = kmeans.clusters)
  return(retList)

}

####--------------------------------------------------------------
## sample.binomial.bmm.component.mean:  Sample (scalar) means in a particular 
##                                      dimension from the posterior of a 
##                                      particular component in the 
##                                      binomial mixture model
##
## Inputs:  
##
## a:  a (scalar) shape parameters for the beta posterior distributions over a
##     particular component's mu parameter (in a particular dimension),
## b:  a (scalar) shape parameters for the beta posterior distributions over a
##     particular component's mu parameter (in a particular dimension)
## num.samples:  the number of samples to return
##
## Outputs:
## 
## a vector of length num.samples holding means sampled from the posterior
## distribution of the component and dimension of the binomial mixture model
## indicated by the input parameters.

sample.binomial.bmm.component.mean <- function(a, b, num.samples) 
{
  # NB:  We could do these integrals analytically, but use this
  # interface for consistency.

  samples.means <- rbeta(num.samples, a, b)

  return(samples.means)

} # End sample.binomial.bmm.component.mean

####--------------------------------------------------------------
## binomial.bmm.narrowest.mean.interval.about.centers:  Return the narrowest (Bayesian
##    credible) interval of the means that include the centers and
##    a fixed proportion of the probability distribution of the means
##    (e.g., .68 or .95).
##
## Inputs:  
##
## a:  a N.c x D matrix holding the _initial_ values of the 
##     shape parameters a for the beta prior distributions over the 
##     mu parameters, i.e., Beta(mu[k,m]; a0[k,m], b0[k,m])  
## b:  a N.c x D matrix holding the _initial_ values of the 
##     shape parameters b for the beta prior distributions over the 
##     mu parameters, i.e., Beta(mu[k,m]; a0[k,m], b0[k,m])  
## proportion:  the percentage of the distribution of each mean to include
##              in the interval (e.g., .68 or .95)
##
## Outputs:
## 
## a list with values lb and ub, where each are D x N.c
## matrices holding the respective lower or upper bounds of the intervals,
## and centers, which is a D x N.c matrix holding the empirical means sampled
## from the posterior distribution

binomial.bmm.narrowest.mean.interval.about.centers <- function(a, b, proportion) 
{

  D <- dim(a)[2]
  N.c <- dim(a)[1]
      
  num.samples <- 1000

  lb <- matrix(data = 0, nrow=D, ncol=N.c)
  ub <- matrix(data = 0, nrow=D, ncol=N.c)
  centers <- matrix(data = 0, nrow=D, ncol=N.c)

  for(i in 1:N.c){
    for(l in 1:D){
      samples.proportions <- sample.binomial.bmm.component.mean(a[i,l], b[i,l], num.samples)

      samples.proportions <- sort(samples.proportions, decreasing=FALSE)
  
      centers[l,i] <- mean(samples.proportions)

      bounds <- define.narrowest.interval.including.pt(samples.proportions, centers[l,i], proportion)
      lb[l,i] <- bounds[1]
      ub[l,i] <- bounds[2]
    }
  }
  
  retList <- list("lb" = lb, "ub" = ub, "centers" = centers)
  return(retList)

} # End binomial.bmm.narrowest.mean.interval.about.centers 

####--------------------------------------------------------------
## binomial.bmm.component.posterior.predictive.density: Calculate the posterior
##   predictive density in a single dimension for a single component.
##
## Inputs:
##
## x:  the integer number of successes at which to evaluate the posterior
##     predictive density
## eta:  the integer total number of trials at which to evaluate the posterior
##     predictive density
## a:  a (scalar) shape parameter of the Beta distribution, i.e.,
##     Beta(x; a, b)
## b:  a (scalar) shape parameter of the Beta distribution, i.e.,
##     Beta(x; a, b)
## pi:  a (scaler) mixing coefficient giving the weight of the component
## Outputs:
##
## the posterior preditive density at x for the specified component

binomial.bmm.component.posterior.predictive.density <- function(x, eta, a, b, pi)
{
  y <- pi * choose(eta, x) * exp(lbeta(x + a, eta - x + b) - lbeta(a, b))
  return(y)
}

####--------------------------------------------------------------
## binomial.bmm.posterior.predictive.density: Calculate the posterior
##   predictive density in a single dimension across all components
##
## Inputs:
##
## x:  the integer number of successes at which to evaluate the posterior
##     predictive density
## eta:  the integer total number of trials at which to evaluate the posterior
##     predictive density
## a:  an N.c vector holding the shape parameter for each of the N.c
##     components
## b:  an N.c vector holding the shape parameter for each of the N.c
##     components
## pi:  an N.c vector holding the mixing coefficients giving the weight 
##      of each of the N.c components
## Outputs:
##
## the posterior preditive density at x summed across all components

binomial.bmm.posterior.predictive.density <- function(x, eta, a, b, pi)
{
  N.c <- length(pi)
  y <- 0
  for(k in 1:N.c) {
    y <- y + binomial.bmm.component.posterior.predictive.density(x, eta, a[k], b[k], pi[k])
  }
  return(y)
}

####--------------------------------------------------------------
## generate.count.data.set:  Generate a data set with clusters of 
##    Binomial-derived counts
##
## Generate N.c clusters of successes, where each dimension is
## independently generated by a binomial distribution.
##
## Inputs:
##
## N.c:  the number of clusters
## N:  the total number of points (which will roughly proportionally
##     allocated across the N.c clusters)
## D:  the number of dimensions of the data set.
##
## Outputs:
##
## a list with slot "successes" holding an N x D matrix of successes (counts)
## a slot "total.trials" holding an N x D matrix of total trials (counts)
##

generate.count.data.set <- function(N.c, N, D)
{
  m <- matrix(data = 0, nrow=N, ncol=D)
  tot.counts <- 100
  tot.count.matrix <- matrix(data = tot.counts, nrow=N, ncol=D)

  num.items.per.cluster <- floor(N/N.c)
  nxt.indx <- 1
  u <- 0
  v <- 0
  for(i in 1:N.c) {
    for(j in 1:D) {
      # Generate a random center for this cluster in this dimension
      u <- runif(1, min=10, max=1000)
      v <- i*runif(1, min=10, max=1000)
      m[nxt.indx:(nxt.indx+num.items.per.cluster-1),j] <- rbinom(num.items.per.cluster, size=tot.counts, prob=rbeta(1, u, v))
    }
    nxt.indx <- nxt.indx + num.items.per.cluster
  }

  # If N does not evenly divide N.c, generate the remaining points
  if(nxt.indx <= N) {
    for(j in 1:D) {
      m[nxt.indx:(nxt.indx+(N-nxt.indx)),j] <- rbinom(N+1-nxt.indx, size=tot.counts, prob=rbeta(1, u, v))
    }
  }

  retList <- list("successes" = m, "total.trials" = tot.count.matrix)
  return(retList)
}

####--------------------------------------------------------------
## binomial.bmm.plot.2d:  Plot data highlighted according to clustering
##
## Plot the data with each item's respective color/symbol indicating
## its cluster.  Also display "1-sigma" contour intervals for the
## standard error of the mean for the standard error.
##
## Inputs:
## successess:  an N x D matrix with rows holding the number of successes
##              in each dimension.
##              All entries are assumed to be counts--not proportions.
## total.trials:  an N x D matrix with rows holding the number of trials
##                in each dimension.
##                All entries are assumed to be counts--not proportions.
## a:  a N.c x D matrix holding the values of the 
##     shape parameters a for the beta posterior distributions over the 
##     mu parameters, i.e., Beta(mu[k,m]; a0[k,m], b0[k,m])  
## b:  a N.c x D matrix holding the values of the 
##     shape parameters b for the beta posterior distributions over the 
##     mu parameters, i.e., Beta(mu[k,m]; a0[k,m], b0[k,m])  
## r:  the N x N.c matrix of responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## title:  plot title
## xlab:  x label
## ylab:  y label
##
## Outputs:
##
## ggplot object displaying the clustered data and "1-sigma" contour intervals
## for the standard error of the mean.

binomial.bmm.plot.2d <- function(successes, total.trials, a, b, r, title, xlab, ylab)
{
  N <- dim(successes)[1]
  N.c <- dim(a)[1]

  X <- successes / total.trials

  # width = 1 std dev
  suppressPackageStartupMessages(library("NORMT3")) # for erf
  width <- as.double(erf(1/sqrt(2)))  
 
  width <- sqrt(width)

  # Calculate standard error of the means
  SEM.res <- binomial.bmm.narrowest.mean.interval.about.centers(a, b, width)
  SEM.centers <- t(SEM.res$centers)
  SEMs.lb <- t(SEM.res$lb)
  SEMs.ub <- t(SEM.res$ub)

  # Make sure none of the clusters = 0 or 1, which would be used for black.
  # Let's reserve that for highlighting.

  clusters <- rep(0, N)
  for(n in 1:N) {
    max.cluster <- 0
    max.assignment <- -1
    for(k in 1:N.c) {
      if ( r[n,k] > max.assignment ) {
        max.assignment <- r[n,k]
        # + 2 ensures none of the clusters has number 0 or 1,
        # which would make it black when cluster is used as colour
        max.cluster <- ( k + 2 )
      }
    }
    clusters[n] <- max.cluster
  }

  proportions <- data.frame(x=X[,1], y=X[,2], row.names=NULL, stringsAsFactors=NULL)

  centers <- data.frame(x=SEM.centers[,1], y=SEM.centers[,2], row.names=NULL, stringsAsFactors=NULL)

  # R CMD check does not understand that the right-hand side of the = in
  # aes commands is referring, not to global or local variables, but to
  # verbatim names in the data argument.  Get around this by creating
  # variables with those names--which is irrelevant to the ggplot code.
  x <- NULL
  y <- NULL
  g <- ggplot(data = proportions, aes(x=x, y=y)) + ggtitle(title) + xlab(xlab) + ylab(ylab) + geom_point(data = proportions, aes(x=x, y=y), shape=clusters, colour=clusters) + geom_point(data = centers, aes(x=x, y=y), size=3, colour='blue') 

  g <- g + theme_bw() + theme(axis.line = theme_segment(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) 
  
  # Add the "1-sigma" contour intervals for the standard error _of the mean_
  for(k in 1:N.c) {
    # i.e., provide the 1-sigma confidence interval of q(mu_k) =
    # \int q(mu_k, lambda_k) dlambda_k.
   
    xc <- SEMs.lb[k,1] + ((SEMs.ub[k,1] - SEMs.lb[k,1])/2)
    yc <- SEMs.lb[k,2] + ((SEMs.ub[k,2] - SEMs.lb[k,2])/2)
  
    ell <- my.ellipse(hlaxa = ((SEMs.ub[k,1] - SEMs.lb[k,1])/2), hlaxb = ((SEMs.ub[k,2] - SEMs.lb[k,2])/2), xc = xc, yc = yc)
  
    df_ell <- cbind(data.frame(x=ell$x, y=ell$y))
  
    g <- g + geom_path(data=df_ell, aes(x=x,y=y), colour="blue", linetype=2)
  }
  
  xmin <- -0.05
  xmax <-  1.05

  g <- g + coord_cartesian(xlim=c(xmin,xmax), ylim=c(xmin,xmax))

  return(g)

} # End binomial.bmm.plot.2d

## End binomial mixture model



## Begin gaussian mixture model

####--------------------------------------------------------------
## gaussian.bmm.fixed.num.components: Use a variational Bayesian approach to fit 
## a mixture of gaussian distributions to proportion data, without dropping
## any components/clusters.  To instead automatically determine the number of
## components, use gaussian.bmm, which invokes this function.
## This implements the derivations described in
##
## Pattern Recognition and Machine Learning
## Christopher M. Bishop
##
## Inputs:
## X:  an N x D matrix with rows being the items to cluster.
##     All entries are assumed to be proportions (i.e., between 0 and 1).
##     Notice that there are no summation restrictions--i.e., proportions
##     do not sum to unity across an item's dimensions.
## N.c:  the number of components/clusters to attempt
## r:  the N x N.c matrix of initial responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## m:  a N.c x D matrix holding the _initial_ centers of the gaussians.
## alpha:  a vector with D components holding the _initial_ values of the
##         parameters of the Dirichlet distribution over the mixing 
##         coefficients pi.  
## beta: a D-vector holding the _initial_ values that scale the precision
## nu:  a D-vector holding the _initial_ values of the degrees of
##      freedom parameter to the Wishart distribution
## W: a list holding the _initial_ D x D symmetric, positive definite matrix
##    parameters to the Wishart distribution, for each of the N.c components
## m0, alpha0, beta0, nu0, W0:  the hyperparameters corresponding to the
##                              above initial values (and with the same
##                              respective matrix/vector dimensionality).
## convergence.threshold:  minimum absolute relative difference between 
##                         variational lower bounds across consecutive
##                         iterations to reach converge.
## max.iterations:  maximum number of iterations to attempt
## verbose:  output progress in terms of mixing coefficient (expected) values
##           if 1.
## Outputs: a list with the following entries
## retVal:  0 indicates successful convergence; -1 indicates a failure
##          to converge.
## m:  a N.c x D matrix holding the _converged_ centers of the gaussians.
## alpha:  a vector with D components holding the _converged_ values of the
##         parameters of the Dirichlet distribution over the mixing 
##         coefficients pi.  
## beta: a D-vector holding the _converged_ values that scale the precision
## nu:  a D-vector holding the _converged_ values of the degrees of
##      freedom parameter to the Wishart distribution
## W: a list holding the _converged_ D x D symmetric, positive definite matrix
##    parameters to the Wishart distribution, for each of the N.c components
## r:  the N x N.c matrix of responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## num.iterations:  the number of iterations required to reach convergence.
## ln.rho:  an N x N.c matrix holding the ln[rho], as defined in eqn (32).
## E.lnpi:  the D-vector holding the values E[ln pi]


gaussian.bmm.fixed.num.components <- function(X, N.c, r, m, alpha, beta, nu, W, m0, alpha0, beta0, nu0, W0, convergence.threshold = 10^-4, max.iterations = 10000, verbose = 0)
{
  N <- dim(X)[1]
  D <- dim(X)[2]

  E.pi.prev <- rep(0, N.c)

  identity.matrix <- matrix(data=0, nrow=D, ncol=D)
  for(d in 1:D) {
    identity.matrix[d,d] <- 1
  }

  # suppressPackageStartupMessages(library("Matrix")) # for lu
  
  W0inv <- W0
  for(k in 1:N.c) {
    # OLD WAY (is best?)
    W0inv[[k]] <- solve(W0[[k]])
    # Take the inverse of W0 by first doing an LU decomposition ...
    #W0.lu <- lu(W0[[k]])
    # ... and using that to define the inverse
    #W0inv[[k]] <- as.matrix(solve(W0.lu))
    delta <- 10^-5    
    if ( any(abs( ( W0[[k]] %*% W0inv[[k]] ) - identity.matrix ) >= delta ) ) {
      cat("W0: \n")
      print(W0[[k]])
      cat("W0inv: \n")
      print(W0inv[[k]])
      stop("Inversion failed")
    }
  }

  # Ensure things don't get too small.  There are sometimes
  # problems with many initial clusters--I assume we divide
  # by r = zero somewhere.
  r <- r + 10^-9
  
  # Bishop 10.51
  Nk <- colSums(r)

  if ( any(Nk <= 0) ) {
    stop("Nk <= 0: Die!\n")
  }

  # Bishop 10.52
  xbar <- matrix(data=0, nrow=N.c, ncol=D)
  for(k in 1:N.c) {
    if ( Nk[k] != 0.0 ) {
      xbar[k,] <- ( 1 / Nk[k] ) * t(X) %*% r[,k]
    }
  }

  # Bishop 10.53
  Sk <- list(length=N.c)
  for(k in 1:N.c) {
    Sk[[k]] <- matrix(data=0, nrow=D, ncol=D)
    if ( Nk[k] != 0.0 ) {
      for(n in 1:N) {
        mean.vec <- X[n,] - xbar[k,]
        Sk[[k]] <- Sk[[k]] + ( r[n,k] * ( mean.vec %o% mean.vec ) )
      }
      Sk[[k]] <- Sk[[k]] / Nk[k]
    }
  }

  iteration <- 0

  # Apply variational bayesian approach to gaussian mixture modeling
  # until convergence

  lb.prev <- -Inf

  while(TRUE) {  

    # M step

    # Bishop eqn 10.58
    alpha <- alpha0 + Nk
    if ( verbose == TRUE ) { cat("alpha = ", alpha, "\n") }
  
    # Bishop eqn 10.60
    beta <- beta0 + Nk
    if ( verbose == TRUE ) { cat("beta = ", beta, "\n") }
  
    # Bishop eqn 10.61
    m <- ( 1 / beta ) * ( beta0 * m0 + Nk * xbar )
    if ( verbose == TRUE ) { cat("m = ", m, "\n") }
  
    # Bishop eqn 10.62
    # NB: solve(W0) == W0^-1 
    W <- W0inv   # Just get the shape of W right.
    for(k in 1:N.c) {
      if ( Nk[k] == 0.0 ) { 
        W[[k]] <- W0[[k]]
        next
      } 
      Winv <- W0inv[[k]] + ( Nk[k] * Sk[[k]] ) + ( ((beta0[k] * Nk[k])/(beta0[k] + Nk[k])) * ( xbar[k,] - m0[k,] ) %o% ( xbar[k,] - m0[k,] ) )
      # Compute inverse by first doing LU decomposition ...
      #Winv.lu <- lu(Winv)
      # ... and then taking inverse based on that decomp.
      #W[[k]] <- as.matrix(solve(Winv.lu))
      W[[k]] <- solve(Winv)
      if ( any(abs( ( W[[k]] %*% Winv ) - identity.matrix ) >= delta ) ) {
        cat("W: \n")
        print(W[[k]])
        cat("Winv: \n")
        print(Winv)
        stop("Inversion failed\n")
      }
    }
  
    # Bishop eqn 10.63
    # Following commented out line is as writtin in Bishop, but I think 
    # this is a typo
    # nu <- nu0 + colSums(r) + 1  
    nu <- nu0 + Nk

    # E step
  
    # Bishop eqn 10.64
    E.quadratic <- matrix(data=0, nrow=N, ncol=N.c)
    for(k in 1:N.c) {
      for(n in 1:N) {
        mean.vec <- X[n,] - m[k,]
        temp <- nu[k] * ( W[[k]] %*% mean.vec )
        val <- mean.vec %*% temp
        val <- val + ( D / beta[k] )
        E.quadratic[n,k] <- val
      }
    }
    if ( verbose == TRUE ) { cat("E.quadratic = ", head(E.quadratic), "\n") }
  
    # Bishop eqn 10.65
    # E.ln.lambda is kappalog in VBmix notation
    E.ln.lambda <- rep(0, N.c)
    gamma.mat <- matrix(data=0, nrow=N.c, ncol=D)
    for(k in 1:N.c) {
      for(d in 1:D) {
        val <- ( nu[k] + 1 - d ) / 2.0
        gamma.mat[k,d] <- digamma(val)
      }
    }
  
    for(k in 1:N.c) {
      val <- sum(gamma.mat[k,])
      val <- val + D * log(2)
      val <- val + determinant(W[[k]], logarithm=TRUE)$modulus
      E.ln.lambda[k] <- val
    }
    if ( verbose == TRUE ) { cat("E.ln.lambda = ", E.ln.lambda, "\n") }
  
    # Bishop eqn 10.66
    # E.ln.pi is pilog in VBmix notation
    E.ln.pi <- digamma(alpha) - digamma(sum(alpha))
    if ( verbose == TRUE ) { cat("E.ln.pi = ", E.ln.pi, "\n") }
  
    # Bishop eqn 10.69: E.pi <- ( alpha0 + Nk ) / ( sum(alpha0) + N )
    # Bishop eqn 10.83
    new.E.pi <- ( 1 / N ) * Nk
  
    # Evaluate responsibilities

    # Bishop eqn 10.46
    pi.matrix <- matrix(data=E.ln.pi, nrow=N, ncol=N.c, byrow=TRUE)
    kappa.matrix <- matrix(data=E.ln.lambda, nrow=N, ncol=N.c, byrow=TRUE)
  
    kappa.matrix <- 0.5 * kappa.matrix
    E.quadratic <- 0.5 * E.quadratic
  
    const.matrix <- matrix(data=(D / 2) * log(2 * pi), nrow=N, ncol=N.c)
  
    ln.rho <- pi.matrix + kappa.matrix - const.matrix - E.quadratic 
      
    if (any(is.na(ln.rho))) {
      print(ln.rho)
      stop("ln.rho is NA")
    }        
  
    # Bishop eqn 10.67
    
    r <- matrix(data = 0, nrow=N, ncol=N.c)
    for(n in 1:N) {
      if(any(is.na(ln.rho[n,]))) {
        r[n,] <- rep(NA, N.c)
        next
      }
      row.sum <- log(sum(exp(ln.rho[n,] - max(ln.rho[n,])))) + max(ln.rho[n,], na.rm=TRUE)
      for(k in 1:N.c) { r[n,k] = exp(ln.rho[n,k] - row.sum) }
    }

    if (any(is.na(r))) {
      print(r)
      stop("r is NA")
    }
    
    if ( verbose == TRUE ) {
      print("    r:")
      print(r)
    }        

    # Ensure things don't get too small.  There are sometimes
    # problems with many initial clusters--I assume we divide
    # by r = zero somewhere.
    r <- r + 10^-9
    
    # Bishop 10.51
    Nk <- colSums(r)
  
    #if(any(Nk <= 0)) { 
    #  cat("Empty cluster:\n")
    #  print(Nk)
    #  # q()
    #}

    if ( verbose == TRUE ) { cat("Nk = ", Nk, "\n") }
  
    # Bishop 10.52
    # xbar is meank in VBmix notation
    xbar <- matrix(data=0, nrow=N.c, ncol=D)
    for(k in 1:N.c) {
      if ( Nk[k] != 0.0 ) {
        xbar[k,] <- ( 1 / Nk[k] ) * t(X) %*% r[,k]
      }
    }
    if ( verbose==TRUE ) { cat("xbar = ", xbar, "\n") }
  
    # Bishop 10.53
    Sk <- list(length=N.c)
    for(k in 1:N.c) {
      Sk[[k]] <- matrix(data=0, nrow=D, ncol=D)
      if ( Nk[k] != 0.0 ) {
        for(n in 1:N) {
          mean.vec <- X[n,] - xbar[k,]
          Sk[[k]] <- Sk[[k]] + ( r[n,k] * ( mean.vec %o% mean.vec ) )
        }
        Sk[[k]] <- Sk[[k]] / Nk[k]
      }
    }
    #if ( verbose == TRUE ) { cat("Sk = ", Sk, "\n") }

    if ( verbose == TRUE) { print(new.E.pi) }
    E.pi <- new.E.pi

    iteration <- iteration + 1
    
    # For efficiency, don't calculate bound on each step
    if ( ( iteration %% 10 ) != 0 ) { next }
    
    # Compute lower bound
    
    # Compute p(X|Z,mu,lambda) Eqn 10.71
    pxlog <- 0
    for(k in 1:N.c) {
      tmp <- E.ln.lambda[k]
      tmp <- tmp - D / beta[k]
      prod <- nu[k] * ( Sk[[k]] %*% W[[k]] )
      tmp <- tmp - sum(diag(prod))
      mean.vec <- xbar[k,] - m[k,]
      temp <- nu[k] * ( W[[k]] %*% mean.vec )
      val <- mean.vec %*% temp
      tmp <- tmp - val
      tmp <- tmp - D * log(2 * pi)
      tmp <- 0.5 * tmp * Nk[k]
      pxlog <- pxlog + tmp
    }
  
    # Compute p(Z|pi) Eqn 10.72
    # NB: this is element-wise multiplication, not matrix-matrix mult
    tab <- r * pi.matrix
    # sum all of the elements of the matrix
    pzlog <- sum(tab)
  
    # Compute p(pi) Eqn 10.73
    ppilog <- lgamma(sum(alpha0)) - sum(lgamma(alpha0)) + sum((alpha0 - 1) * E.ln.pi)
  
    # Compute p(mu, lambda) Eqn 10.74
    pmulog <- 0.5 * D * sum(log(beta0 / (2 * pi)))
    pmulog <- pmulog + 0.5 * sum(E.ln.lambda)
    pmulog <- pmulog - 0.5 * D * sum(beta0 / beta) 
    for(k in 1:N.c) {
      mean.vec <- m[k,] - m0[k,]
      temp.vec <- beta0[k] * nu[k] * ( W[[k]] %*% mean.vec )
      val <- mean.vec %*% temp.vec
      pmulog <- pmulog - val
    }
    # K * ln(B(W0, nu0)) term.  But we allow component-specific
    # W0 and nu0.  So, this becomes an iteration over the K 
    # components (and we lose the K prefactor).
    for(k in 1:N.c) {
      val <- - ( nu0[k] / 2 ) * determinant(W0[[k]], logarithm=TRUE)$modulus
      val <- val - nu0[k] * D * log(2) / 2.0
      val <- val - D * ( D - 1 ) * log(pi) / 4.0
      for(d in 1:D) {
        val <- val - lgamma((nu0[k] + 1 - d) / 2.0)
      }
      pmulog <- pmulog + val
    }
  
    pmulog <- pmulog + 0.5 * sum( (nu0 - D - 1) * E.ln.lambda )
  
    for(k in 1:N.c) {
      prod <- W0inv[[k]] %*% W[[k]]
      pmulog <- pmulog - 0.5 * nu[k] * sum(diag(prod))
    }
  
    # Compute q(Z) Eqn 10.75
    qzlog <- 0
    for(n in 1:N) {
      for(k in 1:N.c) {
        if ( r[n,k] != 0 ) {
          qzlog <- qzlog + r[n,k] * log(r[n,k])
        }
      }
    }
   
    # Compute q(pi) Eqn 10.76
    qpilog <- lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * E.ln.pi)
  
    # Compute q(mu, lambda) Eqn 10.77
    qmulog <- 0.5 * sum(E.ln.lambda) 
    qmulog <- qmulog + 0.5 * D * sum(beta / (2 * pi))
    qmulog <- qmulog - 0.5 * N.c * D
    # Calculate H[q(lambda_k)] term
    for(k in 1:N.c) {
      # Calculate E[ln|lambda_k|]
      L <- 0
      for(d in 1:D) {
        L <- L + digamma( (nu[k] + 1 - d) / 2.0 )
      }
      L <- L + D * log(2)
      W.det <- determinant(W[[k]], logarithm=TRUE)$modulus
      L <- L + W.det
  
      blog <- - ( nu[k] / 2 ) * W.det
      blog <- blog - nu[k] * D * log(2) / 2.0
      blog <- blog - D * ( D - 1 ) * log(pi) / 4.0
      for(d in 1:D) {
        blog <- blog - lgamma((nu[k] + 1 - d) / 2.0)
      }
      val <- - ( (nu[k] - D - 1) / 2.0 ) * L + nu[k] * D / 2.0 - blog
      qmulog <- qmulog - val
    }  
  
    # Calculate the lower bound
    lb <- pxlog + pzlog + ppilog + pmulog - qzlog - qpilog - qmulog

    if ( lb.prev > lb ) {
      cat(sprintf("lb decreased from %f to %f!\n", lb.prev, lb))
      # q(status=-1)
    }
    if ( abs( ( lb - lb.prev ) / lb ) < convergence.threshold ) { break }

    lb.prev <- lb

  } # End inner while(TRUE)

  retList <- list("retVal" = 0, "m" = m, "alpha" = alpha, "beta" = beta, "nu" = nu, "W" = W, "r" = r, "num.iterations" = iteration, "ln.rho" = ln.rho, "E.lnpi" = E.ln.pi, "E.pi" = E.pi)

  return(retList)
} # End gaussian.bmm.fixed.num.components function

####--------------------------------------------------------------
## gaussian.bmm: Use a variational Bayesian approach to fit a mixture of 
## gaussian distributions to count data.
##
## NB: Clustering is first run to convergence using 
## gaussian.bmm.fixed.num.components.
## If any components have probability less than pi.threshold, they are
## discarded and the clustering is run to convergence again.
##
## Inputs:
## X:  an N x D matrix with rows being the items to cluster.
##     All entries are assumed to be proportions (i.e., between 0 and 1).
##     Notice that there are no summation restrictions--i.e., proportions
##     do not sum to unity across an item's dimensions.
## N.c:  the number of components/clusters to attempt
## r:  the N x N.c matrix of initial responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## m:  a N.c x D matrix holding the _initial_ centers of the gaussians.
## alpha:  a vector with D components holding the _initial_ values of the
##         parameters of the Dirichlet distribution over the mixing 
##         coefficients pi.  
## beta: a D-vector holding the _initial_ values that scale the precision
## nu:  a D-vector holding the _initial_ values of the degrees of
##      freedom parameter to the Wishart distribution
## W: a list holding the _initial_ D x D symmetric, positive definite matrix
##    parameters to the Wishart distribution, for each of the N.c components
## m0, alpha0, beta0, nu0, W0:  the hyperparameters corresponding to the
##                              above initial values (and with the same
##                              respective matrix/vector dimensionality).
## convergence.threshold:  minimum absolute relative difference between 
##                         variational lower bounds across consecutive
##                         iterations to reach converge.
## max.iterations:  maximum number of iterations to attempt
## verbose:  output progress in terms of mixing coefficient (expected) values
##           if 1.
## pi.threshold:  discard any cluster with weight/mixing coefficient less
##                than pi.threshold _following_ convergence.  
## Outputs: a list with the following entries
## retVal:  0 indicates successful convergence; -1 indicates a failure
##          to converge.
## m:  a N.c x D matrix holding the _converged_ centers of the gaussians.
## alpha:  a vector with D components holding the _converged_ values of the
##         parameters of the Dirichlet distribution over the mixing 
##         coefficients pi.  
## beta: a D-vector holding the _converged_ values that scale the precision
## nu:  a D-vector holding the _converged_ values of the degrees of
##      freedom parameter to the Wishart distribution
## W: a list holding the _converged_ D x D symmetric, positive definite matrix
##    parameters to the Wishart distribution, for each of the N.c components
## r:  the N x N.c matrix of responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## num.iterations:  the number of iterations required to reach convergence.
## ln.rho:  an N x N.c matrix holding the ln[rho], as defined in eqn (32).
## E.lnpi:  the D-vector holding the values E[ln pi]

gaussian.bmm <- function(X, N.c, r, m, alpha, beta, nu, W, m0, alpha0, beta0, nu0, W0, convergence.threshold = 10^-4, max.iterations = 10000, verbose = 0, pi.threshold = 10^-2)
{
  
  total.iterations <- 0 
  D <- dim(m)[2]

  while(TRUE) {

    if(verbose){
      print(r)
    }

    bmm.res <- gaussian.bmm.fixed.num.components(X, N.c, r, m, alpha, beta, nu, W, m0, alpha0, beta0, nu0, W0, convergence.threshold, max.iterations, verbose)
    
    if(bmm.res$retVal != 0) {
      stop("Failed to converge!\n")
    }

    m <- bmm.res$m
    alpha <- bmm.res$alpha
    beta <- bmm.res$beta
    nu <- bmm.res$nu
    W <- bmm.res$W
    E.pi <- bmm.res$E.pi

    total.iterations <- total.iterations + bmm.res$num.iterations

    ln.rho <- bmm.res$ln.rho
    E.lnpi <- bmm.res$E.lnpi

    do.inner.iteration <- FALSE

    apply.min.items.condition <- TRUE

    if((apply.min.items.condition == TRUE) & (N.c > 1)) {
      non.zero.indices <- E.pi > pi.threshold
      N = length(X[,1])

      if ( any(non.zero.indices==FALSE) ) {

        do.inner.iteration <- TRUE

        numeric.indices <- (1:N.c)

        E.pi <- E.pi[non.zero.indices]
        E.lnpi <- E.lnpi[non.zero.indices]

        N.c <- length(E.pi)
        r <- matrix(r[,non.zero.indices], nrow=N, ncol=N.c)
        ln.rho <- matrix(ln.rho[,non.zero.indices], nrow=N, ncol=N.c)
        # Need to renormalize r--do it gently.
        for(n in 1:N) {
           if(any(is.na(ln.rho[n,]))) {
             r[n,] <- rep(NA, N.c)
             next
           }
          
           row.sum <- log(sum(exp(ln.rho[n,] - max(ln.rho[n,])))) + max(ln.rho[n,])
           for(k in 1:N.c) { r[n,k] = exp(ln.rho[n,k] - row.sum) }
        }

        m <- matrix(m[non.zero.indices,], nrow=N.c, ncol=D)
        alpha <- alpha[non.zero.indices,drop=FALSE]
        beta <- beta[non.zero.indices,drop=FALSE]
        nu <- nu[non.zero.indices,drop=FALSE]
        W <- W[non.zero.indices,drop=FALSE]                        
        m0 <- matrix(m0[non.zero.indices,], nrow=N.c, ncol=D)
        alpha0 <- alpha0[non.zero.indices,drop=FALSE]
        beta0 <- beta0[non.zero.indices,drop=FALSE]
        nu0 <- nu0[non.zero.indices,drop=FALSE]
        W0 <- W0[non.zero.indices,drop=FALSE]                        
      }
    } # End apply.min.items.condition

    if(do.inner.iteration == FALSE) { break }

  }

  retList <- list("retVal" = 0, "m" = m, "alpha" = alpha, "beta" = beta, "nu" = nu, "W" = W, "r" = r, "num.iterations" = total.iterations, "ln.rho" = ln.rho, "E.lnpi" = E.lnpi, "E.pi" = E.pi)

  return(retList)

} # End gaussian.bmm function


####--------------------------------------------------------------
## init.gaussian.bmm.hyperparameters:  Initialize the gaussian mixture model
##                                     hyperparameters (to be passed to 
##                                     gaussian.bmm)
##
## NB:  This provides what should be a generally reasonable initialization of
##      hyperparameters.  However, better results may be obtained by
##      tuning these in an application-specific manner.
##
## Inputs:
## X:  an N x D matrix with rows being the items to cluster.
##     All entries are assumed to be proportions (i.e., between 0 and 1).
##     Notice that there are no summation restrictions--i.e., proportions
##     do not sum to unity across an item's dimensions.
## N.c:  the number of components/clusters to attempt
##
## Outputs:
## m0:  a N.c x D matrix holding the prior centers of the gaussians.
## alpha0:  a vector with D components holding the hyperparameter values of the
##          parameters of the Dirichlet distribution over the mixing 
##          coefficients pi.  
## beta0: a D-vector holding the prior values that scale the precision
## nu0:  a D-vector holding the prior values of the degrees of
##      freedom parameter to the Wishart distribution
## W0: a list holding the prior D x D symmetric, positive definite matrix
##    parameters to the Wishart distribution, for each of the N.c components
init.gaussian.bmm.hyperparameters <- function(X, N.c) 
{
  D <- dim(X)[2]

  # Dirichlet hyperparameter
  # This gives very little prior weight to the existence of each
  # cluster--hence, any cluster that "survives" the iteration will
  # need to be supported by the data.
  alpha0 <- rep(0.001, N.c)
  
  # _Prior_ cluster centers
  m0 <- matrix(data=0, ncol=D, nrow=N.c)
  
  # We have three parameters beta0, nu0, and W0.  We will try to satisfy 3
  # goals:
  # (1) Make the prior of the precision broad 
  # (2) Keep beta0 as nearly as possible to 1--i.e., so that the precision
  #     matrix is coupled to the mu vector and can't go screwy without
  #     affecting the mean.
  # (3) Keep beta0 * E[precision] = beta0 * nu0 * W0 = effective.precision
  # Accomplish these respectively via:
  # (1) Set W0 large
  # (2) Set nu0 as small as possible (as permitted in keeping the
  #     Wishart finite), so that when we guarantee (3) beta0 will be
  #     as large (near 1) as possible.
  # (3) Set beta0 = effective.precision / ( nu0 * W0 )
  
  # Set the desired expected value, effective.precision, of the precision of the 
  # gaussians. This implies an expected value, effective.sigma, of the std deviation,
  # since precision Lambda = Sigma^{-1}.
  # Also, set a desired variance on the sigma.
  effective.precision <- 1000
  #effective.precision <- 10000  
  effective.sigma <- sqrt(1/effective.precision)
  
  delta <- 10^-7
  delta <- 10^-5
  
  # These desired values determine nu0 and W0.
  nu0 <- rep(max(1,D - 1 + delta), N.c)
  #nu0 <- rep(1, N.c)

  W0 <- list(length=N.c)
  for(k in 1:N.c) {
    W0[[k]] <- matrix(data=0, nrow=D, ncol=D)
  }
  
  # NB: we are choosing W0 diagonal and homoscedastic.
  for(d in 1:D) {
    for(k in 1:N.c) {
      W0[[k]][d,d] <- 10^4
    }
  }
  
  beta0 <- rep(0, N.c)
  for(k in 1:N.c) {
    beta0[k] <- effective.precision/(W0[[k]][1,1]*nu0[k])
  }

  retList <- list("m0" = m0, "alpha0" = alpha0, "beta0" = beta0, "nu0" = nu0, "W0" = W0)
  return(retList)

} # End init.gaussian.bmm.hyperparameters


####--------------------------------------------------------------
## init.gaussian.bmm.parameters:  Initialize the gaussian mixture model 
##                                parameters (to be passed to gaussian.bmm)
##
## Initialize parameters such that expected proportions have the values
## determined by an initial k-means clustering.
##
## NB:  This provides what should be a generally reasonable initialization of
##      hyperparameters.  However, better results may be obtained by
##      tuning these in an application-specific manner.
##
##
## Inputs:
## X:  an N x D matrix with rows being the items to cluster.
##     All entries are assumed to be proportions (i.e., between 0 and 1).
##     Notice that there are no summation restrictions--i.e., proportions
##     do not sum to unity across an item's dimensions.
## N.c:  the number of components/clusters to attempt
## m0:  a N.c x D matrix holding the prior centers of the gaussians.
## alpha0:  a vector with D components holding the hyperparameter values of the
##          parameters of the Dirichlet distribution over the mixing 
##          coefficients pi.  
## beta0: a D-vector holding the prior values that scale the precision
## nu0:  a D-vector holding the prior values of the degrees of
##      freedom parameter to the Wishart distribution
## W0: a list holding the prior D x D symmetric, positive definite matrix
##    parameters to the Wishart distribution, for each of the N.c components
##
## Outputs:
## m:  a N.c x D matrix holding the _initial_ centers of the gaussians.
## alpha:  a vector with D components holding the _initial_ values of the
##         parameters of the Dirichlet distribution over the mixing 
##         coefficients pi.  
## beta: a D-vector holding the _initial_ values that scale the precision
## nu:  a D-vector holding the _initial_ values of the degrees of
##      freedom parameter to the Wishart distribution
## W: a list holding the _initial_ D x D symmetric, positive definite matrix
##    parameters to the Wishart distribution, for each of the N.c components
## r:  the N x N.c matrix of initial responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## kmeans.clusters:  an N-vector giving the assignment of each of the N
##                   items to a cluster, as determined by kmeans.
## kmeans.centers:  an N.c x D matrix holding the centers of the N.c
##                  clusters/components determined by kmeans
init.gaussian.bmm.parameters <- function(X, N.c, m0, alpha0, beta0, nu0, W0, verbose=0)
{
  N <- dim(X)[1]
  D <- dim(X)[2]

  # Initialize the responsibilities r_ni with K-means.
  kmeans.out <- kmeans(X, N.c, nstart=1000)
  kmeans.clusters <- kmeans.out$cluster
  kmeans.centers <- kmeans.out$centers

  # Sort the centers to ease comparison across runs.  No we sort the
  # clustered results.
  # Besides, be careful!  This upset the relationship between kmeans.clusters
  # and kmeans.centers!  
  # Sort the centers to ease comparison across runs
  #kmeans.centers <- kmeans.centers[do.call(order, lapply(1:ncol(kmeans.centers), function(i) kmeans.centers[, i])), ]
  if(verbose){
    cat("kmeans initialization:\n")
    write.table(kmeans.centers, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
  }
  r <- matrix(data=0, nrow=N, ncol=N.c)
  for(i in 1:N) {
    r[i,kmeans.clusters[i]]<- 1
  }

  # NB: bmm.gaussian actually computes m, alpha, beta, nu, and W
  # before using them, so we don't need to initialize.  It does
  # use r before computing it, so that does need to be initialized.
  
  # _Initial_ cluster centers
  m <- kmeans.centers

  alpha <- alpha0
  beta <- beta0
  nu <- nu0
  W <- W0
  
  retList <- list("m" = m, "alpha" = alpha, "beta" = beta, "nu" = nu, "W" = W, "r" = r, "kmeans.centers" = kmeans.centers, "kmeans.clusters" = kmeans.clusters)
  return(retList)

}

####--------------------------------------------------------------
## gaussian.bmm.calculate.posterior.predictive.precision: Return the
##    precision of the posterior predictive density.
## Inputs:  
##
## m:  a N.c x D matrix holding the _converged_ centers of the gaussians.
## alpha:  a vector with D components holding the _converged_ values of the
##         parameters of the Dirichlet distribution over the mixing 
##         coefficients pi.  
## beta: a D-vector holding the _converged_ values that scale the precision
## nu:  a D-vector holding the _converged_ values of the degrees of
##      freedom parameter to the Wishart distribution
## W: a list holding the _converged_ D x D symmetric, positive definite matrix
##    parameters to the Wishart distribution, for each of the N.c components
## Outputs:
##
## a list holding the precision matrices of the posterior predictive
## densities for each component.
gaussian.bmm.calculate.posterior.predictive.precision <- function(m, alpha, beta, nu, W)
{
  # Bishop eqn 10.82
  N.c <- dim(m)[1]
  D <- dim(m)[2]

  L <- list(length=N.c)
  for(k in 1:N.c) {
    L[[k]] <- ( ( nu[k] + 1 - D ) * beta[k] ) * W[[k]] / ( 1 + beta[k] )
  }

  return(L)
}

####--------------------------------------------------------------
## gaussian.bmm.narrowest.mean.interval.about.centers:  Return the narrowest (Bayesian
##    credible) interval of the means that include the centers and
##    a fixed proportion of the probability distribution of the means
##    (e.g., .68 or .95).
##
## Inputs:  
##
## m:  a N.c x D matrix holding the _converged_ centers of the gaussians.
## alpha:  a vector with D components holding the _converged_ values of the
##         parameters of the Dirichlet distribution over the mixing 
##         coefficients pi.  
## beta: a D-vector holding the _converged_ values that scale the precision
## nu:  a D-vector holding the _converged_ values of the degrees of
##      freedom parameter to the Wishart distribution
## W: a list holding the _converged_ D x D symmetric, positive definite matrix
##    parameters to the Wishart distribution, for each of the N.c components
## proportion:  the percentage of the distribution of each mean to include
##              in the interval (e.g., .68 or .95)
##
## Outputs:
## 
## a list with values lb and ub, where each are D x N.c
## matrices holding the respective lower or upper bounds of the intervals,
## and centers, which is a D x N.c matrix holding the empirical means sampled
## from the posterior distribution

gaussian.bmm.narrowest.mean.interval.about.centers <- function(m, alpha, beta, nu, W, proportion) 
{
  if(abs(proportion - .68) > .01) {
    stop("gaussian.bmm.narrowest.mean.interval.about.centers only implemented for 1SD = 68%.", proportion, "not supported\n")
  }

  D <- dim(m)[2]
  N.c <- dim(m)[1]

  # Add the "1-sigma" contour intervals for the standard error _of the mean_
  # i.e., provide the 1-sigma confidence interval of q(mu_k) =
  # \int q(mu_k, lambda_k) dlambda_k.
  # As described by Gelman, Carlin, Stern, and Rubin
  # (Bayesian Data Analysis; 2nd edition; page 88), this is
  # also a multivariate student t.  It has degrees of freedom
  # nu[k] + 1 - D (as above), but precision
  # ( ( nu[k] + 1 - D ) * beta[k] ) * W[[k]], i.e.,
  # without the ( 1 + beta[k] ) denominator
  L <- list(length=N.c)
  for(k in 1:N.c) {
    L[[k]] <- ( ( nu[k] + 1 - D ) * beta[k] ) * W[[k]] 
  }

  lb <- matrix(data = 0, nrow=D, ncol=N.c)
  ub <- matrix(data = 0, nrow=D, ncol=N.c)
  centers <- matrix(data = 0, nrow=D, ncol=N.c)

  for(k in 1:N.c) {
    # i.e., provide the 1-sigma confidence interval of the posterior
    # predictive density (Bishop eqn 10.81)--a Student t.
    Lkinv = solve(L[[k]])
    # According to Bishop (B.70), covariance matrix = ( nu / ( nu - 2 ) ) * Lambda^-1
    cov.matrix <- ( ( nu[k] + 1 - D ) / ( nu[k] + 1 - D - 2 ) ) * Lkinv
    for(l in 1:D){
      centers[l,k] <- m[k,l]
      lb[l,k] <- centers[l,k] - cov.matrix[l,l]
      ub[l,k] <- centers[l,k] + cov.matrix[l,l]
    }
  }
  
  retList <- list("lb" = lb, "ub" = ub, "centers" = centers)
  return(retList)

} # End gaussian.bmm.narrowest.mean.interval.about.centers 

####--------------------------------------------------------------
## bmm.narrowest.proportion.interval.about.centers:  Return the narrowest 
##    interval of proportions from the posterior distribution that
##    includes their centers a fixed proportion of the posterior distribution
##    (e.g., .68 or .95).
##
## Inputs:  
## m:  a N.c x D matrix holding the _converged_ centers of the gaussians.
## alpha:  a vector with D components holding the _converged_ values of the
##         parameters of the Dirichlet distribution over the mixing 
##         coefficients pi.  
## beta: a D-vector holding the _converged_ values that scale the precision
## nu:  a D-vector holding the _converged_ values of the degrees of
##      freedom parameter to the Wishart distribution
## W: a list holding the _converged_ D x D symmetric, positive definite matrix
##    parameters to the Wishart distribution, for each of the N.c components
## proportion:  the percentage of the distribution of each mean to include
##              in the interval (e.g., .68 or .95)
##
## Outputs:
## 
## a list with values lb and ub, where each are D x N.c
## matrices holding the respective lower or upper bounds of the intervals,
## and centers, which is a D x N.c matrix holding the empirical means sampled
## from the posterior distribution

gaussian.bmm.narrowest.proportion.interval.about.centers <- function(m, alpha, beta, nu, W, proportion) 
{
  if(abs(proportion - .68) > .01) {
    stop("gaussian.bmm.narrowest.mean.interval.about.centers only implemented for 1SD = 68%.", proportion, "not supported\n")
  }

  L <- gaussian.bmm.calculate.posterior.predictive.precision(m, alpha, beta, nu, W)

  D <- dim(m)[2]
  N.c <- dim(m)[1]

  lb <- matrix(data = 0, nrow=D, ncol=N.c)
  ub <- matrix(data = 0, nrow=D, ncol=N.c)
  centers <- matrix(data = 0, nrow=D, ncol=N.c)

  for(k in 1:N.c) {
    # i.e., provide the 1-sigma confidence interval of the posterior
    # predictive density (Bishop eqn 10.81)--a Student t.
    Lkinv = solve(L[[k]])
    # According to Bishop (B.70), covariance matrix = ( nu / ( nu - 2 ) ) * Lambda^-1    
    cov.matrix <- ( ( nu[k] + 1 - D ) / ( nu[k] + 1 - D - 2 ) ) * Lkinv
    for(l in 1:D){
      centers[l,k] <- m[k,l]
      lb[l,k] <- centers[l,k] - cov.matrix[l,l]
      ub[l,k] <- centers[l,k] + cov.matrix[l,l]
    }
  }
  
  retList <- list("lb" = lb, "ub" = ub, "centers" = centers)
  return(retList)

} # End gaussian.bmm.narrowest.proportion.interval.about.centers 


## Student t distribution
my.dt <- function(x, mu, Lambda, nu)
{
  return( exp(lgamma((nu+1)/2) - lgamma(nu/2)) * sqrt(Lambda/(pi*nu)) * ((1+((Lambda*((x-mu)^2))/(nu)))^(-((nu+1)/2))) )
}


####--------------------------------------------------------------
## gaussian.bmm.component.posterior.predictive.density: Calculate the posterior
##   predictive density in a single dimension for a single component.
##
## Inputs:
##
## x:  a proportion D-vector
## k:  the component at which to evaluate the density
## m:  a N.c x D matrix holding the _converged_ centers of the gaussians.
## alpha:  a vector with D components holding the _converged_ values of the
##         parameters of the Dirichlet distribution over the mixing 
##         coefficients pi.  
## beta: a D-vector holding the _converged_ values that scale the precision
## nu:  a D-vector holding the _converged_ values of the degrees of
##      freedom parameter to the Wishart distribution
## W: a list holding the _converged_ D x D symmetric, positive definite matrix
##    parameters to the Wishart distribution, for each of the N.c components
## rooti:  a list holding the inverses of the cholesky root of the posterior predictive
##    covariance matrix (covariance matrix = inverse of precision matrix)
##    rooti = backsolve(chol(Lkinv), diag(2))
## Outputs:
##
## the posterior preditive density at x for the specified component

gaussian.bmm.component.posterior.predictive.density <- function(x, k, m, alpha, beta, nu, W, L)
{
  D <- dim(m)[2]

  if(D == 1) {
    y <- ( alpha[k] / sum(alpha) ) * my.dt(x, m[k], L[[k]], nu[k] + 1 - D)
    return(y)
  }
  stop("Untested lndMvst")
  #suppressPackageStartupMessages(library("bayesm")) # for lndMvst
  #y <- ( alpha[k] / sum(alpha) ) * lndMvst(x, nu=(nu[k]+1-D), mu=m[k,], rooti=rooti[[k]])
  
  return(y)
}

####--------------------------------------------------------------
## gaussian.bmm.posterior.predictive.density: Calculate the posterior
##   predictive density in a single dimension across all components
##
## Inputs:
##
## x:  a proportion D-vector
## m:  a N.c x D matrix holding the _converged_ centers of the gaussians.
## alpha:  a vector with D components holding the _converged_ values of the
##         parameters of the Dirichlet distribution over the mixing 
##         coefficients pi.  
## beta: a D-vector holding the _converged_ values that scale the precision
## nu:  a D-vector holding the _converged_ values of the degrees of
##      freedom parameter to the Wishart distribution
## W: a list holding the _converged_ D x D symmetric, positive definite matrix
##    parameters to the Wishart distribution, for each of the N.c components
## rooti:  a list holding the inverses of the cholesky root of the posterior predictive
##    covariance matrix (covariance matrix = inverse of precision matrix)
##    rooti = backsolve(chol(Lkinv), diag(2))
## Outputs:
##
## the posterior preditive density at x summed across all components

gaussian.bmm.posterior.predictive.density <- function(x, m, alpha, beta, nu, W, rooti)
{
  y <- 0
  N.c <- length(nu)
  for(k in 1:N.c) {
    y <- y + gaussian.bmm.component.posterior.predictive.density(x, k, m, alpha, beta, nu, W, rooti)
  }
  return(y)
}

####--------------------------------------------------------------
## gaussian.bmm.plot.1d:  Overlay the posterior predictive density on a histogram
##                        of the data.
##
## Inputs:
## X:  an N x D matrix with rows being the items to cluster.
##     All entries are assumed to be proportions (i.e., between 0 and 1).
##     Notice that there are no summation restrictions--i.e., proportions
##     do not sum to unity across an item's dimensions.
## m:  a N.c x D matrix holding the _converged_ centers of the gaussians.
## alpha:  a vector with D components holding the _converged_ values of the
##         parameters of the Dirichlet distribution over the mixing 
##         coefficients pi.  
## beta: a D-vector holding the _converged_ values that scale the precision
## nu:  a D-vector holding the _converged_ values of the degrees of
##      freedom parameter to the Wishart distribution
## W: a list holding the _converged_ D x D symmetric, positive definite matrix
##    parameters to the Wishart distribution, for each of the N.c components
## r:  the N x N.c matrix of responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## title:  plot title
## xlab:  x label
## ylab:  y label
##
## Outputs:
##
## ggplot object displaying the clustered data and "1-sigma" contour intervals
## for the standard error of the mean.

gaussian.bmm.plot.1d <- function(X, m, alpha, beta, nu, W, r, title, xlab, ylab, bin.width = 0.025)
{
  N.c <- dim(m)[1]
  D <- dim(X)[2]

  proportions <- data.frame(x=X[,1], row.names=NULL, stringsAsFactors=NULL)
  
  # Generate (x,y) values of the posterior predictive density
  n <- 1000
  y <- rep.int(0, n)
  # Don't evaluate at x=0 or x=1, which will blow up
  #x <- seq(1/n, 1-(1/n), length=n)
  x <- seq(min(proportions$x), max(proportions$x), length=n)
  ym <- matrix(data=0,nrow=N.c,ncol=n)
  num.iterations <- 1000

  
  rooti <- list(length=N.c)
  L <- gaussian.bmm.calculate.posterior.predictive.precision(m, alpha, beta, nu, W)
  for (k in 1:N.c) {
    Lkinv = solve(L[[k]])
    rooti[[k]] <- backsolve(chol(Lkinv), diag(2))
  }
  
  for (k in 1:N.c) {
    
    for (i in 1:n) {

      # Evaluate posterior probability at x.
      #ym[k,i] <- gaussian.bmm.component.posterior.predictive.density(x, k, m, alpha, beta, nu, W, rooti)
      ym[k,i] <- ( alpha[k] / sum(alpha) ) * my.dt(x[i], m[k], L[[k]], nu[k] + 1 - D)
      #ym[k,i] <- ( alpha[k] / sum(alpha) ) * dnorm(x[i], m[k], sqrt(1/L[[k]]))
  
      y[i] <- y[i] + ym[k,i]
    }
  }

  max.posterior.density <- max(y)
  # Overlay the histogram of the data on to the beta mixture model fit.

  # Set max posterior to max of splinefun.
  limits <- data.frame(x=c(min(x),max(x)))

  f <- splinefun(x,y)
  max.posterior.density <- max(unlist(lapply(seq(from=limits$x[1],to=limits$x[2],by=10^-3), f)))

  g <- ggplot(data = proportions, aes(x)) + ggtitle(title) + xlab(xlab) + ylab(ylab)

  g <- g + theme_bw() + theme(axis.line = theme_segment(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) 

  num.breaks <- ceiling(1/bin.width) + 1
  breaks <- unlist(lapply(0:(num.breaks-1), function(x) x/(num.breaks-1)))
  #g <- g + geom_histogram(data=proportions, mapping=aes(x), fill="white", colour="black", breaks=breaks)
  g <- g + geom_histogram(data=proportions, mapping=aes(x), fill="white", colour="black")

  tmp <- print(g)
  max.density <- max(tmp[["data"]][[1]]$ymax)

  scale <- (max.density/max.posterior.density)

  vertical.offset <- .1 * max.posterior.density

  for(k in 1:N.c) {
      
    f <- splinefun(x,scale*ym[k,])

    g <- g + stat_function(data = limits, fun=f, mapping=aes(x), lty="dashed")

  }

  f <- splinefun(x,scale*y)

  g <- g + stat_function(data = limits, fun=f, mapping=aes(x))
  
  #xmin <- -0.05
  #xmax <-  1.05
  xmin <- min(proportions$x)
  xmax <- max(proportions$x)

  g <- g + coord_cartesian(ylim=c(0, max.density*1.1), xlim=c(xmin,xmax))

  return(g)
} # End gaussian.bmm.plot.1d


####--------------------------------------------------------------
## gaussian.bmm.plot.2d:  Plot data highlighted according to clustering
##
## Plot the data with each item's respective color/symbol indicating
## its cluster.  Also display "1-sigma" contour intervals for the
## standard error of the mean for the standard error.
##
## Inputs:
## X:  an N x D matrix with rows being the items to cluster.
##     All entries are assumed to be proportions (i.e., between 0 and 1).
##     Notice that there are no summation restrictions--i.e., proportions
##     do not sum to unity across an item's dimensions.
## m:  a N.c x D matrix holding the _converged_ centers of the gaussians.
## alpha:  a vector with D components holding the _converged_ values of the
##         parameters of the Dirichlet distribution over the mixing 
##         coefficients pi.  
## beta: a D-vector holding the _converged_ values that scale the precision
## nu:  a D-vector holding the _converged_ values of the degrees of
##      freedom parameter to the Wishart distribution
## W: a list holding the _converged_ D x D symmetric, positive definite matrix
##    parameters to the Wishart distribution, for each of the N.c components
## r:  the N x N.c matrix of responsibilities, with r[n, nc] giving the
##     probability that item n belongs to component nc
## title:  plot title
## xlab:  x label
## ylab:  y label
##
## Outputs:
##
## ggplot object displaying the clustered data and "1-sigma" contour intervals
## for the standard error of the mean.

gaussian.bmm.plot.2d <- function(X, m, alpha, beta, nu, W, r, title, xlab, ylab)
{
  N <- dim(X)[1]
  N.c <- dim(m)[1]

  print(N.c)
  
  # width = 1 std dev
  suppressPackageStartupMessages(library("NORMT3")) # for erf
  width <- as.double(erf(1/sqrt(2)))  
 
  # width <- sqrt(width)

  # Calculate standard error of the means
  SEM.res <- gaussian.bmm.narrowest.mean.interval.about.centers(m, alpha, beta, nu, W, width) 
  SEM.centers <- t(SEM.res$centers)
  SEMs.lb <- t(SEM.res$lb)
  SEMs.ub <- t(SEM.res$ub)

  # Make sure none of the clusters = 0 or 1, which would be used for black.
  # Let's reserve that for highlighting.

  clusters <- rep(0, N)
  for(n in 1:N) {
    max.cluster <- 0
    max.assignment <- -1
    for(k in 1:N.c) {
      if ( r[n,k] > max.assignment ) {
        max.assignment <- r[n,k]
        # + 2 ensures none of the clusters has number 0 or 1,
        # which would make it black when cluster is used as colour
        max.cluster <- ( k + 2 )
      }
    }
    clusters[n] <- max.cluster
  }

  proportions <- data.frame(x=X[,1], y=X[,2], row.names=NULL, stringsAsFactors=NULL)

  centers <- data.frame(x=SEM.centers[,1], y=SEM.centers[,2], row.names=NULL, stringsAsFactors=NULL)

  
  # R CMD check does not understand that the right-hand side of the = in
  # aes commands is referring, not to global or local variables, but to
  # verbatim names in the data argument.  Get around this by creating
  # variables with those names--which is irrelevant to the ggplot code.
  x <- NULL
  y <- NULL
  g <- ggplot(data = proportions, aes(x=x, y=y)) + ggtitle(title) + xlab(xlab) + ylab(ylab) + geom_point(data = proportions, aes(x=x, y=y), shape=clusters, colour=clusters) + geom_point(data = centers, aes(x=x, y=y), size=3, colour='blue') 

  g <- g + theme_bw() + theme(axis.line = theme_segment(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) 
  
  # Add the "1-sigma" contour intervals for the standard error _of the mean_
  for(k in 1:N.c) {
    # i.e., provide the 1-sigma confidence interval of q(mu_k) =
    # \int q(mu_k, lambda_k) dlambda_k.
   
    xc <- SEMs.lb[k,1] + ((SEMs.ub[k,1] - SEMs.lb[k,1])/2)
    yc <- SEMs.lb[k,2] + ((SEMs.ub[k,2] - SEMs.lb[k,2])/2)
  
    ell <- my.ellipse(hlaxa = ((SEMs.ub[k,1] - SEMs.lb[k,1])/2), hlaxb = ((SEMs.ub[k,2] - SEMs.lb[k,2])/2), xc = xc, yc = yc)
  
    df_ell <- cbind(data.frame(x=ell$x, y=ell$y))
  
    g <- g + geom_path(data=df_ell, aes(x=x,y=y), colour="blue", linetype=2)
  }
  
  xmin <- -0.05
  xmax <-  1.05

  g <- g + coord_cartesian(xlim=c(xmin,xmax), ylim=c(xmin,xmax))

  return(g)

} # End gaussian.bmm.plot.2d

## End gaussian mixture model



