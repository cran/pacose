trk <-
function(ind,x) {
  ni <- length(ind)
  n <- nrow(x)
  if (ni==1) {
    return(-(1/((n-1)*var(x[,ind+1])))**2)
  }else{
    mattemp <- solve((n-1)*cov(x[,ind+1]))
    return(-1/2*tr(mattemp**2) - 1/2*(tr(mattemp))**2 )
  }
}

