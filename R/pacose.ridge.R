pacose.ridge <-
function (X, gg, lambda = NULL, scale = FALSE, k = 10, verbose = FALSE, cv.method="CV") {
    if (is.null(lambda) == TRUE) {
        ss <- seq(-10, -1, length = 1000)
        ss <- 10^ss
        n <- nrow(X)
        nn <- n - floor(n/k)
        lambda <- ss * nn * ncol(X)
    }
    n <- nrow(X)
    p <- ncol(X)
    pv <- pvar.shrink(X,verbose=F)
    X <- scale(X, scale = scale)
    A <- get.adjacency(gg)
    B <- matrix(0, nrow = p, ncol = p)
    lambda.opt <- rep(0, p)
    #pv <- rep(0, p)
    if (verbose) cat(paste("Performing local ridge regressions\n"))
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
              if (cv.method == "CV" ) lambda.opt.i <- ridge.cv(Xi, yi, lambda = lambda, scale = scale, plot.it = FALSE, k = k)$lambda.opt
              if (cv.method == "HKB") lambda.opt.i <- lm.ridge(yi~Xi, lambda = 0, scale = scale)$kHKB
              
              rr <- lm.ridge(yi ~ Xi, scale = scale, lambda = lambda.opt.i)
              
              B[i, noti ] <- coef(rr)[-1]
              lambda.opt[i] <- lambda.opt.i
              #pv[i] <- var(yi-(Xi%*%coef(rr)[-1]+coef(rr)[1]))
          } else {
              lmo <- lm(yi~Xi)
              B[i, noti ] <- lmo$coefficients[2]
              lambda.opt[i] <- 0
              #pv[i] <- var(lmo$residuals)
          }
        }  
    }
    if (verbose) cat("\n")
    # print(lambda.opt)
    pcor <- Beta2parcor(B, verbose = verbose)
    invcov <- pcor2invcov(pcor,pv)
    return(list(pcor = pcor,invcov=invcov))
}

