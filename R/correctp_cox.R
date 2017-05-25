
correctp.cox=
function (x, y, eta, K, kappa, select, fit, verbose = FALSE)
{
    force(K)
    if (min(eta) < 0 | max(eta) >= 1) {
        if (max(eta) == 1) {
            stop("eta should be strictly less than 1!")
        }
        if (length(eta) == 1) {
            stop("eta should be between 0 and 1!")
        }
        else {
            stop("eta should be between 0 and 1! \n  Choose appropriate range of eta!")
        }
    }
    if (max(K) > ncol(x)) {
        stop("K cannot exceed the number of predictors! Pick up smaller K!")
    }
    if (max(K) >= nrow(x)) {
        stop("K cannot exceed the sample size! Pick up smaller K!")
    }
    if (min(K) <= 0 | !all(K%%1 == 0)) {
        if (length(K) == 1) {
            stop("K should be a positive integer!")
        }
        else {
            stop("K should be a positive integer! \n  Choose appropriate range of K!")
        }
    }
    if (kappa > 0.5 | kappa < 0) {
        if (verbose) {
            cat("kappa should be between 0 and 0.5! kappa=0.5 is used. \n\n")
        }
        kappa <- 0.5
    }
    if (select != "pls2" & select != "simpls") {
        if (verbose) {
            cat("Invalid PLS algorithm for variable selection.\n")
        }
        if (verbose) {
            cat("pls2 algorithm is used. \n\n")
        }
        select <- "pls2"
    }
    fits <- c("regression", "canonical", "invariant", "classic")
    if (!any(fit == fits)) {
        if (verbose) {
            cat("Invalid PLS algorithm for model fitting\n")
        }
        if (verbose) {
            cat("regression algorithm is used. \n\n")
        }
        fit <- "regression"
    }
    list(K = K, eta = eta, kappa = kappa, select = select, fit = fit)
}

