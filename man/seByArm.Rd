\name{seByArm}
\alias{seByArm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  seByArm: split a SummarizedExperiment by chromosomal arm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
seByArm(SE, GRorGRL = NULL, build = "hg19")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SE}{
%%     ~~Describe \code{SE} here~~
}
  \item{GRorGRL}{
%%     ~~Describe \code{GRorGRL} here~~
}
  \item{build}{
%%     ~~Describe \code{build} here~~
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
function (SE, GRorGRL = NULL, build = "hg19") 
{
    stopifnot(is(SE, "SummarizedExperiment"))
    if (!is.null(GRorGRL)) {
        stopifnot(unique(genome(GRorGRL)) == build)
        stopifnot(class(GRorGRL) \%in\% c("GRanges", "GRangesList"))
        stopifnot(unique(genome(rowData(SE))) == unique(genome(GRorGRL)))
    }
    else {
        data(hg19.by.arm)
        GRorGRL <- hg19.by.arm
    }
    if (is(GRorGRL, "GRanges")) 
        GRL <- split(GRorGRL)
    else GRL <- GRorGRL
    lapply(GRL, function(x) SE[names(subsetByOverlaps(rowData(SE), 
        x)), ])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
