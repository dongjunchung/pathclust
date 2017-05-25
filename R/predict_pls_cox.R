
predict.pls.cox<-
function (object, newdata, scale.X = TRUE, scale.Y = TRUE, ...)
{
    if (missing(newdata))
        stop("No new data available.")
    X <- object$X
    Y <- object$Y
    q <- ncol(Y)
    p <- ncol(X)
    if (length(dim(newdata)) == 2) {
        if (ncol(newdata) != p)
            stop("'newdata' must be a numeric matrix with ncol = ",
                p, " or a vector of length = ", p, ".")
    }
    if (length(dim(newdata)) == 0) {
        if (length(newdata) != p)
            stop("'newdata' must be a numeric matrix with ncol = ",
                p, " or a vector of length = ", p, ".")
        dim(newdata) <- c(1, p)
    }
    ncomp <- object$ncomp
    a <- object$loadings$X
    b <- object$loadings$Y
    c <- object$mat.c
    if (scale.X) {
        means.X <- attr(X, "scaled:center")
    }
    if (scale.Y) {
        means.Y <- attr(Y, "scaled:center")
    }
    if (scale.X) {
        sigma.X <- attr(X, "scaled:scale")
    }
    if (scale.Y) {
        sigma.Y <- attr(Y, "scaled:scale")
    }
    newdata <- as.matrix(newdata)
    ones <- matrix(rep(1, nrow(newdata)), ncol = 1)
    B.hat <- array(0, dim = c(p, q, ncomp))
    Y.hat <- array(0, dim = c(nrow(newdata), q, ncomp))
    t.pred <- array(0, dim = c(nrow(newdata), ncomp))
    for (h in 1:ncomp) {
        W <- a[, 1:h] %*% solve(t(c[, 1:h]) %*% a[, 1:h])
        B <- W %*% drop(t(b[, 1:h]))
        if (scale.Y) {
            B <- scale(B, center = FALSE, scale = 1/sigma.Y)
        }
        if (scale.X) {
            B <- as.matrix(scale(t(B), center = FALSE, scale = sigma.X))
        }
        if (!scale.X) {
            B <- as.matrix(t(B))
        }
        if (scale.X & scale.Y) {
            intercept <- -scale(B, center = FALSE, scale = 1/means.X)
            intercept <- matrix(apply(intercept, 1, sum) + means.Y,
                nrow = 1)
            Y.hat[, , h] <- newdata %*% t(B) + ones %*% intercept
        }
        if (scale.X & !scale.Y) {
            intercept <- -scale(B, center = FALSE, scale = 1/means.X)
            intercept <- matrix(apply(intercept, 1, sum), nrow = 1)
            Y.hat[, , h] <- newdata %*% t(B) + ones %*% intercept
        }
        if (!scale.X & scale.Y) {
            intercept <- -B
            intercept <- matrix(apply(intercept, 1, sum) + means.Y,
                nrow = 1)
            Y.hat[, , h] <- newdata %*% t(B) + ones %*% intercept
        }
        if (!scale.X & !scale.Y) {
            Y.hat[, , h] <- newdata %*% t(B)
        }
        if (!scale.X) {
            t.pred[, h] <- newdata %*% W[, h]
        }
        if (scale.X) {
            t.pred[, h] <- scale(newdata, center = means.X, scale = sigma.X) %*%
                W[, h]
        }
        B.hat[, , h] <- B
    }
    rownames(t.pred) <- rownames(newdata)
    colnames(t.pred) <- paste("dim", c(1:ncomp), sep = " ")
    rownames(Y.hat) <- rownames(newdata)
    colnames(Y.hat) <- colnames(Y)
    return(invisible(list(predict = Y.hat, variates = t.pred, w=W,
        B.hat = B.hat)))
}

