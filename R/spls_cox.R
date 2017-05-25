
spls.cox <- function (x, y, K, eta, kappa = 0.5, select = "pls2", fit = "regression",
    scale.x = TRUE, scale.y = FALSE, eps = 1e-04, maxstep = 100,
    trace = FALSE)
{
    force(K)
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    ip <- c(1:p)
    y <- as.matrix(y)
    q <- ncol(y)
    one <- matrix(1, 1, n)
    mu <- one %*% y/n
    y <- scale(y, drop(mu), FALSE)
    meanx <- drop(one %*% x)/n
    x <- scale(x, meanx, FALSE)
    if (scale.x) {
        normx <- sqrt(drop(one %*% (x^2))/(n - 1))
        if (any(normx < .Machine$double.eps)) {
            stop("Some of the columns of the predictor matrix have zero variance.")
        }
        x <- scale(x, FALSE, normx)
    }
    else {
        normx <- rep(1, p)
    }
    if (scale.y) {
        normy <- sqrt(drop(one %*% (y^2))/(n - 1))
        if (any(normy < .Machine$double.eps)) {
            stop("Some of the columns of the response matrix have zero variance.")
        }
        y <- scale(y, FALSE, normy)
    }
    else {
        normy <- rep(1, q)
    }
    betahat <- matrix(0, p, q)
    betamat <- list()
    x1 <- x
    y1 <- y
    type <- correctp.cox(x, y, eta, K, kappa, select, fit)
    eta <- type$eta
    K <- type$K
    kappa <- type$kappa
    select <- type$select
    fit <- type$fit
    if (is.null(colnames(x))) {
        xnames <- c(1:p)
    }
    else {
        xnames <- colnames(x)
    }
    new2As <- list()
    if (trace) {
        cat("The variables that join the set of selected variables at each step:\n")
    }

    for (k in 1:K) {
        Z <- t(x1) %*% y1
        what <- spls.dv(Z, eta, kappa, eps, maxstep)
        A <- unique(ip[what != 0 | betahat[, 1] != 0])
        new2A <- ip[what != 0 & betahat[, 1] == 0]
        xA <- x[, A, drop = FALSE]
        plsfit <- pls.cox(X = xA, Y = y, ncomp = min(k, length(A)),
            mode = fit, scale.X = FALSE, scale.Y = FALSE)

        predplsfit <- predict.pls.cox(plsfit, newdata = xA, scale.X = FALSE,
            scale.Y = FALSE)
        betahat <- matrix(0, p, q)
        betahat[A, ] <- matrix(predplsfit$B.hat[, , plsfit$ncomp],
            length(A), q)
        betamat[[k]] <- betahat
        if (select == "pls2") {
            y1 <- y - predplsfit$predict[, , plsfit$ncomp]
        }
        new2As[[k]] <- new2A
        if (trace) {
            if (length(new2A) <= 10) {
                cat(paste("- ", k, "th step (K=", k, "):\n",
                  sep = ""))
                cat(xnames[new2A])
                cat("\n")
            }
            else {
                cat(paste("- ", k, "th step (K=", k, "):\n",
                  sep = ""))
                nlines <- ceiling(length(new2A)/10)
                for (i in 0:(nlines - 2)) {
                  cat(xnames[new2A[(10 * i + 1):(10 * (i + 1))]])
                  cat("\n")
                }
                cat(xnames[new2A[(10 * (nlines - 1) + 1):length(new2A)]])
                cat("\n")
            }
        }
    }
    if (!is.null(colnames(x))) {
        rownames(betahat) <- colnames(x)
    }
    if (q > 1 & !is.null(colnames(y))) {
        colnames(betahat) <- colnames(y)
    }
    object <- list(x = x, y = y, betahat = betahat, A = A, betamat = betamat,
        new2As = new2As, mu = mu, meanx = meanx, normx = normx,
        normy = normy, eta = eta, K = K, kappa = kappa, select = select,
        fit = fit, projection = NA, plsmod = plsfit, pred=predplsfit)
    class(object) <- "spls"
    object
}

