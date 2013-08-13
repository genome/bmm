\name{bmm.posterior.predictive.density}

\alias{bmm.posterior.predictive.density}

\title{
  bmm.posterior.predictive.density: Calculate the posterior
    predictive density.

}

\description{
  Calculate the posterior predictive density at a point x in a 
  single dimension across all components
}

\usage{
  bmm.posterior.predictive.density(x, mu, alpha, nu, beta, pi, num.samples = 1000)
}

\arguments{
   \item{x}{ 
      the scalar at which to evaluate the posterior predictive density
   }
   \item{mu}{ 
      an N.c vector holding the shape parameter for each of the N.c 
      components, which governs the gamma prior distribution over that 
      component's u parameter.
   }
   \item{alpha}{ 
      an N.c vector holding the rate (i.e., inverse scale) parameter 
      for each of the N.c components, which governs the gamma prior 
      distribution over that component's u parameter.
   }
   \item{nu}{ 
      an N.c vector holding the shape parameter for each of the N.c 
      components, which governs the gamma prior distribution over that 
      component's v parameter.
   }
   \item{beta}{ 
      an N.c vector holding the rate (i.e., inverse scale) parameter 
      for each of the N.c components, which governs the gamma prior 
      distribution over that component's v parameter.
   }
   \item{pi}{ 
      an N.c vector holding the mixing coefficients giving the weight 
      of each of the N.c components
   }
   \item{num.samples}{  
      the number of samples to use in performing the numerical
      evaluation of the predictive density
   }
}

\value{
  the posterior preditive density at x summed across all components
}



