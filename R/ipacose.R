ipacose <-
function(x=x,pc=pc,method="pacose.ridge",cutoff,gr=NULL,k=2, cv.method="CV",adaptive=FALSE) {
  nonstop <- TRUE
  p <- ncol(x)
  iter <- 0
  gr1 <- gr2 <- graph.adjacency((abs(pc) >= cutoff) - diag(p),mode="undirected")
  while(nonstop) {
   iter <- iter+1
   res2 <- do.call(method,list(X=x,gg=gr2,k=k, cv.method=cv.method))
   if (method == "pacose.adalasso") {
      if (adaptive==TRUE){ pcor <- res2$pcor.lasso
      }else{ pcor <- res2$pcor.adalasso }
   }else{
      pcor <- res2$pcor
   }
   #if (is.null(cutoff) & iter==1) {
   #  cutoff <- fdrtool(pc[upper.tri(pc)], statistic = "correlation", plot = F,verbose=F)$param[1]
   #  print("Estimated cutoff = " ,cutoff,"\n")
   #}
   # print(round(c(cutoff,nbri(gr,gr2)/nbr(gr),nbri(gr,gr2)/nbr(gr2)),3))
   if (cutoff==0) {
      grnew <- gr2
   }else{
      grnew <- graph.adjacency( (abs(pcor) >= cutoff) - diag(p) ,mode="undirected")
   }
   nonstop <- nbri(grnew,gr2) != nbru(grnew,gr2)
   gr2 <- grnew
  }
  if (is.null(gr)) {
    sens <- ppv <- NULL
  } else {
    sens <- nbri(gr,gr2)/nbr(gr); ppv <- nbri(gr,gr2)/nbr(gr2)
  }
  return(list(Niter=iter,gr_it=gr2,pcor_it=pcor,sens=sens,ppv=ppv))
}

