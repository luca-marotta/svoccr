\name{getCoOccMatrices}
\alias{getCoOccMatrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
getCoOccMatrices
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCoOccMatrices(sym.comb, co.occ)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sym.comb}{
%%     ~~Describe \code{sym.comb} here~~
}
  \item{co.occ}{
%%     ~~Describe \code{co.occ} here~~
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
function (sym.comb, co.occ) 
{
    mat.strat <- strsplit(sym.comb, split = ",")
    mat.strat <- lapply(mat.strat, strtoi)
    mat.list <- lapply(seq_along(mat.strat), function(i) {
        m <- matrix(0, 3, 3)
        m[as.matrix(co.occ)[mat.strat[[i]], , drop = F]] <- 1
        m
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
