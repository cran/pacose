isComplete <-
function (G, set){
 if (length(set)<=1) {
   Res <- TRUE
 } else {
   B       <- as(G, "matrix")
   B       <- B[set,set]
   diag(B) <- 1
   Res     <- all(B != 0)
 }
 return(Res)
}
