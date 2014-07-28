\name{writeCombMat}
\alias{writeCombMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
writeCombMat
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
writeCombMat(sym.comb, pairs, s = 3, outfile = "./sym_combinations.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sym.comb}{
%%     ~~Describe \code{sym.comb} here~~
}
  \item{pairs}{
%%     ~~Describe \code{pairs} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
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
function (sym.comb, pairs, s = 3, outfile = "./sym_combinations.txt") 
{
    write("keys for the co-occurrences of strategies:\n", file = outfile)
    co.occ <- expand.grid(a1 = seq_len(s), a2 = seq_len(s))
    states.key <- data.frame(co.occ, key = seq_len(nrow(co.occ)))
    suppressWarnings(write.table(states.key, file = outfile, 
        quote = F, row.names = F, sep = "\t", append = T))
    mat <- getCoOccMatrices(names(sym.comb), co.occ)
    write("Strategy Combinations. Please refer to the key up above for the co-occurence code.\n", 
        append = T, file = outfile)
    for (i in seq(sym.comb)) {
        cat("Combination", names(sym.comb)[i], "frequency:", 
            sym.comb[i], "\n", file = outfile, append = T, sep = " ")
        cat("Combination Matrix representation:\n", file = outfile, 
            append = T, sep = " ")
        write.table(mat[[i]], file = outfile, col.names = F, 
            row.names = F, quote = F, append = T)
        cat("Combination edgelist:\n", file = outfile, append = T)
        suppressWarnings(write.table(pairs[[i]], file = outfile, 
            col.names = T, row.names = F, quote = F, sep = "\t", 
            append = T))
        cat("\n-------------------------------------\n", file = outfile, 
            append = T)
    }
    return(0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
