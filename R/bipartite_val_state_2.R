bipartite_val<-function(infile="./input.txt", s=2, outfile="./validation.txt"){
  edgel <- read.table(infile, header=F, sep="\t");
  colnames(edgel) <- c("agent", "day", "start", "end", "activity");
  edgel <- edgel[order(edgel$agent),];
  #### sparse adjacency matrix
  adj <-matrix(0, nrow=length(unique(edgel[,1])), ncol=max(edgel$end), 
               dimnames=list(unique(edgel[,1]), 1:max(edgel$end)));
  #### filling the matrix
  adj[as.matrix(edgel[,1:2])] <- edgel$activity;
  
  #if(!is.matrix(adj)) return("Please note adj must be the bipartite adjacency matrix");
  
  state <- sort(unique(adj[which(adj!=0)]));
  if(!length(state)==s) stop("the number of states in adj does no match input states");
  
  co.occ <- expand.grid(state, state);
  ns <- nrow(co.occ);
  
  #### dataframe with start and end date for every agent
  df <- edgel[!duplicated(edgel[,1]), c("agent", "start", "end")];
  #require(igraph);
  ############### functiont to calculate the p-value after formula 1 in PLOS ###############
  hyper_pvalue <- function(sizeA, sizeB, inter.size, totsize){
    if(inter.size==0) return(1);
    return(1 - phyper(inter.size-1, sizeA, totsize-sizeA, sizeB));
  }
  
  ###### setting names for rows and colums of adj if not present
  #if(is.null(dimnames(adj)[[1]])) rownames(adj)<-paste("A", seq(1:nrow(adj)), sep="");
  #if(is.null(dimnames(adj)[[2]])) colnames(adj)<-paste("B", seq(1:ncol(adj)), sep="");  
  
  ############### constructing the related graph 
  g <- graph.incidence(adj, directed=F, weighted="weight");   
  
  ######### projecting the bipartite graph
  g.proj<-bipartite.projection(g)[[1]];
  
  ###### setting the Bonferroni threshold: The number of tests B1 is equal to the number of pairs 
  ###### of vertices in the projected network. The univariate threshold is 1%.  
  B1 <- vcount(g.proj)*(vcount(g.proj)-1)*0.5*nrow(co.occ);
  B1.thresh <- 0.01/B1;
  
  ###########validation function
  #### a1, a2: id of nodes
  #### sa1, sa2: state of nodes
  #validate_state_pair <- function(a1, a2, sa1, sa2) {
  validate_state_pair <- function(a1, a2) {
    common_adj <- adj[c(a1, a2), max(df[a1,]$start, df[a2,]$start):min(df[a1,]$end, df[a2,]$end), drop=F];
    if(ncol(common_adj)==0) return(1);
    hp <-function(sa1, sa2){
      hyper_pvalue(sum(common_adj[1,]==sa1), sum(common_adj[2,]==sa2), 
                   sum(common_adj[1,]==sa1 & common_adj[2,]==sa2),
                   ncol(common_adj));
    }
    
    pval <- mapply(hp, co.occ[,1], co.occ[,2]);
    m <- matrix(c(rep(a1, ns), rep(a2, ns), pval, 1:ns), ns, 4);
    res <- m[m[,3]<B1.thresh, ,drop=F];
    if(nrow(res)>=1) write.table(res, outfile, col.names=F, row.names=F, quote=F, sep="\t", append=T);    
    return(m);
  } 
  
  #### cycling through all pairs to validate every state
  state.pairs <- t(combn(1:vcount(g.proj), 2));
    
  system.time({el <- mapply(validate_state_pair, state.pairs[,1], state.pairs[,2])});
  
  el <-as.data.frame(matrix(t(el), nrow(co.occ)*ncol(el), 4)); colnames(el) <- c("a1", "a2", "pvalue", "co.occ");
  
  #### creating validated network
  g.val <- graph.data.frame(el[el$pvalue < B1.thresh, c(1,2,4)], directed=F);
  E(g.val)$color <- palette(rainbow(nrow(co.occ)))[E(g.val)$co.occ];
  E(g.val)$weight <- 1; ### to be used for community detection
  ###################################
  #adding color column to co.occ matrix
  co.occ <- data.frame(co.occ, co.occ=1:nrow(co.occ), color=palette(rainbow(nrow(co.occ))));
  
  return(list(co.occ=co.occ, network=g.val, edgelist=el, bonf.thresh=B1.thresh));
}