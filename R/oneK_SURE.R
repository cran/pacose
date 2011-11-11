oneK_SURE <-
function(ind,x,d) {
  n <- nrow(x) ; p <- ncol(x)
  mat1 <- matrix(0,p,p)
  ni <- length(ind)
  if (ni==1) {
    mat1[ind+1,ind+1] <- (n-ni-1-d)/(n-1)*1/(var(x[,ind+1]))
  }else{
    mat1[ind+1,ind+1] <- (n-ni-1-d)*solve((n-1)*cov(x[,ind+1]))
  }
  mat1
}

