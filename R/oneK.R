oneK <-
function(ind,x) {
  n <- nrow(x) ; p <- ncol(x)
  mat1 <- matrix(0,p,p)
  ni <- length(ind)
  if (ni==1) {
    mat1[ind+1,ind+1] <- 1/((n-1)*var(x[,ind+1]))
  }else{
    mat1[ind+1,ind+1] <- solve((n-1)*cov(x[,ind+1]))
  }
  mat1
}

