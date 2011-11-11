omegaMVUE <-
function(x,gr) {
  if (!(is.chordal(gr,fillin=TRUE)$chordal)) {warning("Non decomposable graph");return(matrix(NaN,ncol(x),ncol(x)))}
  m <- get.adjacency(gr,sparse=F) ; colnames(m) <- rownames(m) <- 1:ncol(m)
  jt <- ug.to.jtree(m)
  l1 <- lapply(jt@cliques,function(u) as.numeric(u$vset))
  l2 <- lapply(jt@separators,function(u) as.numeric(u$separator))
  p <- ncol(x)
  KhatMVUE <- Reduce("+",lapply(l1,oneK_MVUE,x)) - Reduce("+",lapply(l2,oneK_MVUE,x))
  return(KhatMVUE=KhatMVUE)
}
