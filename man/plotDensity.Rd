\name{plotDensity}
\alias{plotDensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  plotDensity
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotDensity(x, minSep = 0.1, xlab = "Value", more = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{minSep}{
%%     ~~Describe \code{minSep} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{more}{
%%     ~~Describe \code{more} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, minSep = 0.1, xlab = "Value", more = "", ...) 
{
    x <- x[!is.na(x)]
    res = clusterAndCollapse(x, minSep)
    means = unique(res$values)
    plot(density(x, adjust = 0.5), xlab = xlab, lty = 1, lwd = 2, 
        ..., main = paste0(res$retain, " components (dashed) at minSep=", 
            minSep, more))
    if (res$retain > 0) {
        for (i in seq_along(means)) abline(v = means[i], col = i + 
            2, lty = 3, lwd = 2)
        dens = density(res$values, adjust = 0.5)
        lines(dens$x, dens$y * (1.2/length(means)), lty = 3, 
            lwd = 2)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
