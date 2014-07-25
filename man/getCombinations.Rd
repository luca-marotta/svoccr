\name{getCombinations}
\alias{getCombinations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
getCombinations
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCombinations(edgel, co.occ)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edgel}{
%%     ~~Describe \code{edgel} here~~
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
function (edgel, co.occ) 
{
    require(plyr)
    el.val1 <- plyr::ddply(as.data.frame(edgel), .variables = c("a1", 
        "a2"), summarize, comb = paste(co.occ, collapse = ","))
    sym.co.occ <- match(data.frame(t(co.occ)), data.frame(t(co.occ)[2:1, 
        ]))
    edgel <- data.frame(edgel, sym.co.occ = sym.co.occ[edgel$co.occ])
    el2 <- data.frame(lapply(data.frame(t(edgel[, c("co.occ", 
        "sym.co.occ")])), sort))
    el.val2 <- data.frame(edgel[, 1:2], t(el2))
    colnames(el.val2) <- c("a1", "a2", "co.occ", "sym.co.occ")
    el.val2 <- data.frame(el.val1, sym.comb = plyr::ddply(as.data.frame(el.val2), 
        .variables = c("a1", "a2"), summarize, comb = paste(co.occ, 
            sym.co.occ, collapse = ","))$comb)
    freq <- table(el.val2$sym.comb)
    pairs <- lapply(names(freq), function(i) el.val2[el.val2$sym.comb == 
        i, 1:3])
    names(freq) <- el.val2$comb[match(names(freq), el.val2$sym.comb)]
    return(list(sym.comb = freq, pairs = pairs, co.occ = co.occ))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
