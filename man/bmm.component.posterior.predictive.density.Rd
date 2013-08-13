\name{bmm.component.posterior.predictive.density}

\alias{bmm.component.posterior.predictive.density}

\title{
  bmm.component.posterior.predictive.density: Calculate a component's posterior
     predictive density.
}

\description{
  Calculate the posterior predictive density at x in a single dimension for a 
  single component.
}

\usage{
  bmm.component.posterior.predictive.density(x, mu, alpha, nu, beta, pi, 
                                             num.samples = 1000)
}

\arguments{
   \item{x}{
     the scalar at which to evaluate the posterior predictive density
   }
   \item{mu}{
     a (scalar) shape parameter for the gamma prior distribution over
     a particular component's u parameter.
   }
   \item{alpha}{
     a (scalar) rate (i.e., inverse scale) parameter for the gamma prior 
     distribution over a particular component's u parameter.
   }
   \item{nu}{
     a (scalar) shape parameter for the gamma prior distribution over
     a particular component's v parameter.
   }
   \item{beta}{
     a (scalar) rate (i.e., inverse scale) parameter for the gamma prior 
     distribution over a particular component's v parameter.
   }
   \item{pi}{
     a (scaler) mixing coefficient giving the weight of the component
   }
   \item{num.samples}{
     the number of samples to use in performing the numerical
     evaluation of the predictive density
   }
}

\value{
  the posterior preditive density at x for the specified component
}




