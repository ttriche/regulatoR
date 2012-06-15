\name{keepFeature}
\alias{keepFeature}
\title{
  keepFeature: logical vector of mixture outcomes at a given separation
}
\description{
  Tells you whether or not to keep a feature at a specified separation
}
\usage{
keepFeature(x, minSep = 0.2, G = c(1:5), ...)
}
\arguments{
  \item{x}{
    A matrix
  }
  \item{minSep}{
    An integer (minimum group mean separation)
  }
  \item{G}{
    Numbers of mixture components to try fitting
  }
  \item{\dots}{
    Anything else I forgot 
  }
}
\value{
  A logical vector
}
\author{
  Tim Triche, Jr.
}

\examples{
  # keepFeature()
}
