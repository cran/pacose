omegaSURE <-
function(x,gr) {
  if (!isDecomposable(igraph.to.graphNEL(gr))) {warning("Non decomposable graph");return(matrix(NaN,ncol(x),ncol(x)))}
  struct <- rip(igraph.to.graphNEL(gr))
  l1 <- lapply(struct$cliques,as.numeric)
  l2 <- lapply(struct$separators,function(vec) ifelse(length(vec)==0,NA,as.numeric(vec)) )
  l2 <- l2[!sapply(l2,function(xx) any(is.na(xx)))]
  p <- ncol(x)
  normD2 <- sum((Reduce("+",lapply(l1,oneK,x)) - Reduce("+",lapply(l2,oneK,x)))**2)
  num <-  -2*( sum(sapply(l1,trk,x)) - sum(sapply(l2,trk,x)) )
  d <- num/normD2
  KhatSURE <- Reduce("+",lapply(l1,oneK_SURE,x,d)) - Reduce("+",lapply(l2,oneK_SURE,x,d))
  return(KhatSURE=KhatSURE)
}

