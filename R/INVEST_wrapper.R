INVEST_wrapper <-
function(x,ind,delta=10e-10,itermax = 2000) {
 tempmat <- matrix(10,ncol(x),ncol(x))
 maxiter <- itermax
 indices <- ind
 icovx <- pseudoinverse(cov(x))
 p <- ncol(x)
 invcovx_vecform <- .C("wermuthC", p = as.integer(p), ninteract=as.integer(nrow(indices))
                                  , delta = as.double(delta), error = as.double(1)
                                  , iter = as.integer(0), maxiter = as.integer(maxiter)
                                  , ind = as.integer(t(indices)), icovx = as.double(icovx)
                                  , result = as.double(tempmat),PACKAGE="pacose")$result
 invcovx <-  matrix(invcovx_vecform,p,p)
 error   <- max(abs(invcovx[indices+1]))
 return(list(invcovx=invcovx,error=error))
}

