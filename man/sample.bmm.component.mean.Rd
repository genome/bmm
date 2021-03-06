\name{sample.bmm.component.mean}

\alias{sample.bmm.component.mean}

\title{
  sample.bmm.component.mean:  Sample (scalar) means from a bmm component
}

\description{
  sample.bmm.component.mean:  Sample (scalar) means in a 
                              particular dimension from the 
                              posterior of a particular component in 
                              the Beta mixture model
}

\usage{
  sample.bmm.component.mean(mu, alpha, nu, beta, num.samples)
}

\arguments{
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
  \item{num.samples}{
     the number of samples to return
  }
}

\value{
   a vector of length num.samples holding means sampled from the  
   posterior distribution of the component and dimension of the Beta 
   mixture model indicated by the input parameters.
}



