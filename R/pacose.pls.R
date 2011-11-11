pacose.pls <-
function (X, gg, scale = TRUE, k = 10, Ncomp = NULL, verbose = FALSE, cv.method="CV") {
    n <- nrow(X)
    p <- ncol(X)
    pv <- pvar.shrink(X,verbose=F)
    k <- max(1,floor(k))
    if (k > n) {
        cat(paste("k exceeds the number of observations. Leave-one-out is applied.\n"))
        k <- n
    }
    A <- get.adjacency(gg)
    B <- matrix(0, nrow = p, ncol = p)
    m <- vector(length = p)
    if (verbose) cat(paste("Performing local pls regressions\n"))
    kernel <- FALSE
    if (n < (p - 1)) {
        kernel <- TRUE
    }
    if (verbose) cat(paste("Vertex no "))
    for (i in 1:p) {
        if ((i/10) == floor(i/10)) {
            if (verbose) cat(paste(i, "..."))
        }
        noti <- (1:p)[-which(A[i,] == 0)]
        yi <- X[, i]
        if(length(noti)==0) { #pv[i] <- var(yi)
        } else {
          Xi <- matrix(X[, noti],ncol=length(noti))
          Xi <- X[, noti]
          if (!is.null(dim(Xi))) {
              if (is.null(Ncomp)) {
                     ncomp <- min(n - 1, ncol(Xi))
              } else {
                     ncomp <- Ncomp
              }
              fit <- penalized.pls.cv(Xi, yi, scale = scale, k = k, ncomp = ncomp)
              B[i, noti] <- fit$coefficients
              m[i] <- fit$ncomp.opt
          } else {
              lmo <- lm(yi~Xi)
              B[i, noti ] <- lmo$coefficients[2]
              m[i] <- 0
              #pv[i] <- var(lmo$residuals)
          }
        }
    }
    if (verbose) cat("\n")
    pcor <- Beta2parcor(B, verbose = verbose)
    invcov <- pcor2invcov(pcor,pv)
    return(list(pcor = pcor,invcov=invcov, m = m))
}

