\name{bipartite_val}
\alias{bipartite_val}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
bipartite_val
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bipartite_val(infile = "~/input.txt", s = 2, outfile = "~/validation.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{
%%     ~~Describe \code{infile} here~~
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
function (infile = "~/input.txt", s = 2, outfile = "~/validation.txt") 
{
    edgel <- read.table(infile, header = F, sep = "\t")
    colnames(edgel) <- c("agent", "day", "start", "end", "activity")
    edgel <- edgel[order(edgel$agent), ]
    adj <- matrix(0, nrow = length(unique(edgel[, 1])), ncol = max(edgel$end), 
        dimnames = list(unique(edgel[, 1]), 1:max(edgel$end)))
    adj[as.matrix(edgel[, 1:2])] <- edgel$activity
    state <- sort(unique(adj[which(adj != 0)]))
    if (!length(state) == s) 
        stop("the number of states in adj does no match input states")
    co.occ <- expand.grid(state, state)
    ns <- nrow(co.occ)
    df <- edgel[!duplicated(edgel[, 1]), c("agent", "start", 
        "end")]
    require(igraph)
    hyper_pvalue <- function(sizeA, sizeB, inter.size, totsize) {
        if (inter.size == 0) 
            return(1)
        return(1 - phyper(inter.size - 1, sizeA, totsize - sizeA, 
            sizeB))
    }
    g <- graph.incidence(adj, directed = F, weighted = "weight")
    g.proj <- bipartite.projection(g)[[1]]
    B1 <- vcount(g.proj) * (vcount(g.proj) - 1) * 0.5 * nrow(co.occ)
    B1.thresh <- 0.01/B1
    validate_state_pair <- function(a1, a2) {
        common_adj <- adj[c(a1, a2), max(df[a1, ]$start, df[a2, 
            ]$start):min(df[a1, ]$end, df[a2, ]$end), drop = F]
        if (ncol(common_adj) == 0) 
            return(1)
        hp <- function(sa1, sa2) {
            hyper_pvalue(sum(common_adj[1, ] == sa1), sum(common_adj[2, 
                ] == sa2), sum(common_adj[1, ] == sa1 & common_adj[2, 
                ] == sa2), ncol(common_adj))
        }
        pval <- mapply(hp, co.occ[, 1], co.occ[, 2])
        m <- matrix(c(rep(a1, ns), rep(a2, ns), pval, 1:ns), 
            ns, 4)
        res <- m[m[, 3] < B1.thresh, , drop = F]
        if (nrow(res) >= 1) 
            write.table(res, outfile, col.names = F, row.names = F, 
                quote = F, sep = "\t", append = T)
        return(m)
    }
    state.pairs <- t(combn(1:vcount(g.proj), 2))
    system.time({
        el <- mapply(validate_state_pair, state.pairs[, 1], state.pairs[, 
            2])
    })
    el <- as.data.frame(matrix(t(el), nrow(co.occ) * ncol(el), 
        4))
    colnames(el) <- c("a1", "a2", "pvalue", "co.occ")
    g.val <- graph.data.frame(el[el$pvalue < B1.thresh, c(1, 
        2, 4)], directed = F)
    E(g.val)$color <- palette(rainbow(nrow(co.occ)))[E(g.val)$co.occ]
    E(g.val)$weight <- 1
    co.occ <- data.frame(co.occ, co.occ = 1:nrow(co.occ), color = palette(rainbow(nrow(co.occ))))
    return(list(co.occ = co.occ, network = g.val, edgelist = el, 
        bonf.thresh = B1.thresh))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
