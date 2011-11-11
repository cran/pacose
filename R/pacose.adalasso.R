pacose.adalasso <-
function (X, gg, k = 10, use.Gram = FALSE, both = TRUE, verbose = FALSE, cv.method="CV")  {
    p <- ncol(X)
    X <- scale(X)
    colnames(X) <- 1:p
    A <- get.adjacency(gg)
    B.lasso <- B.adalasso <- matrix(0, nrow = p, ncol = p)
    colnames(B.lasso) <- colnames(B.adalasso) <- 1:p
    pcor.adalasso <- NULL
    pv <- pvar.shrink(X,verbose=F)
    if (verbose) cat(paste("Performing local (adaptive) lasso regressions\n"))
    if (verbose) cat(paste("Vertex no "))
    for (i in 1:p) {
        if ((i/10) == floor(i/10)) {
            if (verbose) cat(paste(i, "..."))
        }
        noti <- (1:p)[-which(A[i,] == 0)]
        #noti <- (1:p)[-i]
        yi <- X[, i]
        Xi <- X[, noti]
        if(length(noti)==0) { #pv[i] <- var(yi)
        } else {
          if (!is.null(dim(Xi))) {
            dummy <- adalasso(Xi, yi, k = k, use.Gram = use.Gram,both = both)
            coefi.lasso <- dummy$coefficients.lasso
            B.lasso[i, noti] <- coefi.lasso
            if (both == TRUE) {
               coefi.adalasso <- dummy$coefficients.adalasso
               B.adalasso[i, noti] <- coefi.adalasso
            }
          } else {
              lmo <- lm(yi~Xi)
              B.lasso[i, noti ] <- lmo$coefficients[2]
              if (both == TRUE) {
                 B.adalasso[i, noti] <- lmo$coefficients[2]
              }
              #pv[i] <- var(lmo$residuals)
          }

        }
    }
    pcor.lasso <- Beta2parcor(B.lasso, verbose = verbose)
    invcov.lasso <- pcor2invcov(pcor.lasso,pv)
    if (both == TRUE) {
        pcor.adalasso <- Beta2parcor(B.adalasso, verbose = verbose)
        invcov.adalasso <- pcor2invcov(pcor.adalasso,pv)
    }
    if (verbose) cat(paste("\n"))
    return(list(pcor.lasso = pcor.lasso, pcor.adalasso = pcor.adalasso,
                invcov.lasso = invcov.lasso, invcov.adalasso = invcov.adalasso))
}

