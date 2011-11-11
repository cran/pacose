INVEST_wrapper <-
function(x,ind,delta=10e-10,itermax = 2000) {
 temp <- matrix(10,ncol(x),ncol(x))
 maxiter <- itermax
 indices <- ind
 icovx <- pseudoinverse(cov(x))
 p <- ncol(x)
 invcovx <-  matrix(.C("INVEST_wrapper", p = as.integer(p), ninteract=as.integer(nrow(indices))
                                  , delta = as.double(delta), error = as.double(1)
                                  , iter = as.integer(0), maxiter = as.integer(maxiter)
                                  , ind = as.integer(t(indices)), icovx = as.double(t(icovx))
                                  , result = as.double(temp),PACKAGE="pacose")$result,p,p
                                  )

 error <- max(abs(invcovx[indices+1]))
 return(list(invcovx=invcovx,error=error))
}

