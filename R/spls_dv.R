#' @import stats

spls.dv=
function (Z, eta, kappa, eps, maxstep)
{
    p <- nrow(Z)
    q <- ncol(Z)
    Znorm1 <- median(abs(Z))
    Z <- Z/Znorm1
    if (q == 1) {
        c <- ust(Z, eta)
    }
    if (q > 1) {
        M <- Z %*% t(Z)
        dis <- 10
        i <- 1
        if (kappa == 0.5) {
            c <- matrix(10, p, 1)
            c.old <- c
            while (dis > eps & i <= maxstep) {
                mcsvd <- svd(M %*% c)
                a <- mcsvd$u %*% t(mcsvd$v)
                c <- ust(M %*% a, eta)
                dis <- max(abs(c - c.old))
                c.old <- c
                i <- i + 1
            }
        }
        if (kappa > 0 & kappa < 0.5) {
            kappa2 <- (1 - kappa)/(1 - 2 * kappa)
            c <- matrix(10, p, 1)
            c.old <- c
            h <- function(lambda) {
                alpha <- solve(M + lambda * diag(p)) %*% M %*%
                  c
                obj <- t(alpha) %*% alpha - 1/kappa2^2
                return(obj)
            }
            if (h(eps) * h(1e+30) > 0) {
                while (h(eps) <= 1e+05) {
                  M <- 2 * M
                  c <- 2 * c
                }
            }
            while (dis > eps & i <= maxstep) {
                if (h(eps) * h(1e+30) > 0) {
                  while (h(eps) <= 1e+05) {
                    M <- 2 * M
                    c <- 2 * c
                  }
                }
                lambdas <- uniroot(h, c(eps, 1e+30))$root
                a <- kappa2 * solve(M + lambdas * diag(p)) %*%
                  M %*% c
                c <- ust(M %*% a, eta)
                dis <- max(abs(c - c.old))
                c.old <- c
                i <- i + 1
            }
        }
    }
    return(c)
}

