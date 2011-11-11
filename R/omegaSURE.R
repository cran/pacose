omegaSURE <-
function(x,gr) {
  if (!(is.chordal(gr,fillin=TRUE)$chordal)) {warning("Non decomposable graph");return(matrix(NaN,ncol(x),ncol(x)))}
  m <- get.adjacency(gr,sparse=F) ; colnames(m) <- rownames(m) <- 1:ncol(m)
  jt <- ug.to.jtree(m)
  l1 <- lapply(jt@cliques,function(u) as.numeric(u$vset))
  l2 <- lapply(jt@separators,function(u) as.numeric(u$separator))
  p <- ncol(x)
  normD2 <- sum((Reduce("+",lapply(l1,oneK,x)) - Reduce("+",lapply(l2,oneK,x)))**2)
  num <-  -2*( sum(sapply(l1,trk,x)) - sum(sapply(l2,trk,x)) )
  d <- num/normD2
  KhatSURE <- Reduce("+",lapply(l1,oneK_SURE,x,d)) - Reduce("+",lapply(l2,oneK_SURE,x,d))
  return(KhatSURE=KhatSURE)
}

