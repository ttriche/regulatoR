\name{clusterAndCollapse}
\alias{clusterAndCollapse}
\title{
  clusterAndCollapse: fit a mixture model with minimum mean separation
}
\description{
  This calls Mclust and then greedily collapses not-well-separated components.
}
\usage{
  clusterAndCollapse(x, minSep, G = c(1:5))
}
\arguments{
  \item{x}{
    A vector or matrix
  }
  \item{minSep}{
    Minimum separation required between means of adjacent groups
  }
  \item{G}{
    Numbers of groups to fit (default is 1 through 5)
  }
}
\details{
  Not yet.
}
\value{
  A list with two elements: 
  $retain is the number of surviving components, 
  $values is the reassigned values for each observation 
}
\author{
  Tim Triche, Jr..
}
\examples{
  
  foo <- c(0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9)
  clusterAndCollapse(foo, minSep=0.4)

}
