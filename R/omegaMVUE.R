omegaMVUE <-
function(x,gr) {
  if (!isDecomposable(igraph.to.graphNEL(gr))) {warning("Non decomposable graph");return(matrix(NaN,ncol(x),ncol(x)))}
  struct <- rip(igraph.to.graphNEL(gr))
  l1 <- lapply(struct$cliques,as.numeric)
  l2 <- lapply(struct$separators,function(vec) ifelse(length(vec)==0,NA,as.numeric(vec)) )
  l2 <- l2[!sapply(l2,function(xx) any(is.na(xx)))]
  p <- ncol(x)
  KhatMVUE <- Reduce("+",lapply(l1,oneK_MVUE,x)) - Reduce("+",lapply(l2,oneK_MVUE,x))
  return(KhatMVUE=KhatMVUE)
}

