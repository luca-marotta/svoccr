################ example script for bipartite_val_state
getCombinations <- function(edgel, co.occ){
  #require(plyr);
  el.val1 <- plyr::ddply(as.data.frame(edgel), .variables=c("a1", "a2"), summarize, comb=paste(co.occ, collapse=","));
  sym.co.occ <- match(data.frame(t(co.occ)), data.frame(t(co.occ)[2:1, ]));
  edgel <- data.frame(edgel, sym.co.occ=sym.co.occ[edgel$co.occ]);
  
  el2 <- data.frame(lapply(data.frame(t(edgel[,c("co.occ", "sym.co.occ")])), sort));
  
  el.val2 <- data.frame(edgel[,1:2], t(el2));
  colnames(el.val2) <- c("a1", "a2","co.occ", "sym.co.occ");
  el.val2 <- data.frame(el.val1, sym.comb=plyr::ddply(as.data.frame(el.val2), .variables=c("a1", "a2"), summarize, comb=paste(co.occ, sym.co.occ, collapse=","))$comb);
  #### results
  ## vector containing unique symmetric combinations of co-occurences and their frequences
  freq <- table(el.val2$sym.comb); 
  ### list with validated pairs and their co-occurrences for every combination
  pairs <- lapply(names(freq), function(i) el.val2[el.val2$sym.comb==i, 1:3]);
  names(freq) <- el.val2$comb[match(names(freq), el.val2$sym.comb)];
  
  return(list(sym.comb=freq, pairs=pairs, co.occ=co.occ)); 
  #return(list(sym.comb=sym.comb, pairs=pairs, co.occ=co.occ)); 
}

getCoOccMatrices <- function(sym.comb, co.occ){
  #### obtaining numeric indices from the combinations
  mat.strat <- strsplit(sym.comb, split=",");
  mat.strat <- lapply(mat.strat, strtoi);
  #use indices to subset the co.occ matrix
  mat.list <- lapply(seq_along(mat.strat), function(i){
    m <- matrix(0, 3, 3);
    m[as.matrix(co.occ)[mat.strat[[i]],, drop=F]] <- 1; 
    m;});   
}

writeCombMat <- function(sym.comb, pairs, s=3, outfile="./sym_combinations.txt"){
  write("keys for the co-occurrences of strategies:\n", file = outfile);
  co.occ <- expand.grid(a1=seq_len(s), a2=seq_len(s))
  states.key <- data.frame(co.occ, key=seq_len(nrow(co.occ)));
  suppressWarnings(write.table(states.key, file=outfile, quote=F, row.names=F, sep="\t", append=T));
  mat <- getCoOccMatrices(names(sym.comb), co.occ);
  write("Strategy Combinations. Please refer to the key up above for the co-occurence code.\n", append=T, file=outfile);    
  for(i in seq(sym.comb)){
    cat("Combination", names(sym.comb)[i], "frequency:", sym.comb[i], "\n", file=outfile, append=T, sep=" ");
    cat("Combination Matrix representation:\n", file=outfile, append=T, sep=" ");
    write.table(mat[[i]], file=outfile, col.names=F, row.names=F, quote=F, append=T);
    cat("Combination edgelist:\n", file=outfile, append=T);
    suppressWarnings(write.table(pairs[[i]], file=outfile, col.names=T, row.names=F, quote=F, sep="\t", append=T));
    cat("\n-------------------------------------\n", file=outfile, append=T);
  }
  
}