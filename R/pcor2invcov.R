pcor2invcov <-
function(A,pv) {
   # if (any(is.nan(A))) return(NA)
   A <- -A ; diag(A) <- -diag(A)
   vec <- 1/pv
   sqrt(diag(vec))%*%A%*%sqrt(diag(vec))
}

