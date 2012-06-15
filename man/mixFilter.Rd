\name{mixFilter}
\alias{mixFilter}
\title{
  mixFilter()
}
\description{
  run clusterAndCollapse but don't blow up if it throws an error
}
\usage{
  mixFilter(x, minSep = 0.2, G = c(1:5), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A vector of observations
  }
  \item{minSep}{
    Minimum separation for mixture group means
  }
  \item{G}{
    Numbers of mixture components to fit
  }
  \item{\dots}{
    Stuff I forgot
  }
}
\value{
  A list with elements $retain and $values, as for clusterAndCollapse
}
\author{
  Tim Triche, Jr.
}
\examples{
  # mixFilter(x, 0.2)
}
