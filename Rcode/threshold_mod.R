threshold.wd.mod <- function (wd, levels = 3:(nlevelsWT(wd) - 1), type = "soft", 
    policy = "sure", by.level = FALSE, value = 0, dev = madmad, 
    boundary = FALSE, verbose = FALSE, return.threshold = FALSE, 
    force.sure = FALSE, cvtol = 0.01, cvmaxits = 500, Q = 0.05, 
    OP1alpha = 0.05, alpha = 0.5, beta = 1, C1 = NA, C2 = NA, 
    C1.start = 100, al.check = TRUE, ...) 
{
    if (verbose == TRUE) 
        cat("threshold.wd:\n")
    if (IsEarly(wd)) {
        ConvertMessage()
        stop()
    }
    if (verbose == TRUE) 
        cat("Argument checking\n")
    ctmp <- class(wd)
    if (is.null(ctmp)) 
        stop("wd has no class")
    else if (ctmp != "wd") 
        stop("wd is not of class wd")
    if (policy != "universal" && policy != "manual" && policy != 
        "probability" && policy != "sure" && policy != "mannum" && 
        policy != "cv" && policy != "fdr" && policy != "op1" && 
        policy != "op2" && policy != "LSuniversal" && policy != 
        "BayesThresh") 
        stop("Only policys are universal, BayesThresh, manual, mannum, sure, LSuniversal, cv, op1, op2 and probability at present")
    if (type != "hard" && type != "soft") 
        stop("Only hard or soft thresholding at  present")
    r <- range(levels)
    if (r[1] < 0) 
        stop("levels out of range, level too small. Minimum level is 0")
    if (r[2] > nlevelsWT(wd) - 1) 
        stop(paste("levels out of range, level too big. Maximum level is", 
            nlevelsWT(wd) - 1))
    if (r[1] > nlevelsWT(wd) - 1) {
        warning("no thresholding done")
        return(wd)
    }
    if (r[2] < 0) {
        warning("no thresholding done")
        return(wd)
    }
    if (al.check == TRUE) 
        if (all(sort(levels) == levels) == FALSE) 
            warning("Entries in levels vector are not ascending. Please check this is what you intend. If so, you can turn this warning off with al.check argument")
    d <- NULL
    n <- 2^nlevelsWT(wd)
    nthresh <- length(levels)
    if (is.complex(wd$D)) {
        stop("Please use cthresh package for complex-valued wavelet shrinkage")
    }
    if (policy == "universal") {
        if (verbose == TRUE) 
            cat("Universal policy...")
        if (by.level == FALSE) {
            if (verbose == TRUE) 
                cat("All levels at once\n")
            for (i in 1:nthresh) d <- c(d, accessD(wd, level = levels[i], 
                boundary = boundary))
            dd <- accessD(wd, level = levels[nthresh], boundary = boundary)
            noise.level <- sqrt(dev(dd))
            nd <- length(d)
            thresh <- sqrt(2 * log(nd)) * noise.level
            if (verbose == TRUE) 
                cat("Global threshold is: ", thresh, "\n")
            thresh <- rep(thresh, length = nthresh)
        }
        else {
            if (verbose == TRUE) 
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            for (i in 1:nthresh) {
                d <- accessD(wd, level = levels[i], boundary = boundary)
                noise.level <- sqrt(dev(d))
                nd <- length(d)
                thresh[i] <- sqrt(2 * log(nd)) * noise.level
                if (verbose == TRUE) 
                  cat("Threshold for level: ", levels[i], " is ", 
                    thresh[i], "\n")
            }
        }
    }
    else if (policy == "LSuniversal") {
        if (verbose == TRUE) 
            cat("Local spectral universal policy...")
        if (by.level == FALSE) {
            if (verbose == TRUE) 
                cat("All levels at once\n")
            for (i in 1:nthresh) d <- c(d, accessD(wd, level = levels[i], 
                boundary = boundary))
            noise.level <- sqrt(dev(d))
            nd <- length(d)
            thresh <- log(nd) * noise.level
            if (verbose == TRUE) 
                cat("Global threshold is: ", thresh, "\n")
            thresh <- rep(thresh, length = nthresh)
        }
        else {
            if (verbose == TRUE) 
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            for (i in 1:nthresh) {
                d <- accessD(wd, level = levels[i], boundary = boundary)
                noise.level <- sqrt(dev(d))
                nd <- length(d)
                thresh[i] <- log(nd) * noise.level
                if (verbose == TRUE) 
                  cat("Threshold for level: ", levels[i], " is ", 
                    thresh[i], "\n")
            }
        }
    }
    else if (policy == "sure") {
        if (type == "hard") 
            stop("Can only do soft thresholding with sure policy")
        if (by.level == FALSE) {
            if (verbose == TRUE) 
                cat("All levels at once\n")
            for (i in 1:nthresh) d <- c(d, accessD(wd, level = levels[i], 
                boundary = boundary))
            dd <- accessD(wd, level = levels[nthresh], boundary = boundary)
            noise.level <- sqrt(dev(dd))
            nd <- length(d)
            neta.d <- (log(nd, base = 2)^(3/2))
            sd2 <- (sum((d/noise.level)^2 - 1)/nd)
            if (verbose == TRUE) {
                cat("neta.d is ", neta.d, "\nsd2 is ", sd2, "\n")
                cat("nd is ", nd, "\n")
                cat("noise.level ", noise.level, "\n")
            }
            if (force.sure == TRUE || sd2 > neta.d/sqrt(nd)) {
                if (verbose == TRUE) {
                  cat("SURE: Using SURE\n")
                }
                thresh <- sure(d/noise.level)
            }
            else {
                if (verbose == TRUE) 
                  cat("SURE: (sparse) using sqrt 2log n\n")
                thresh <- sqrt(2 * log(nd))
            }
            thresh <- rep(thresh * noise.level, length = nthresh)
            if (verbose == TRUE) 
                cat("Global threshold is ", thresh, "\n")
        }
        else {
            if (verbose == TRUE) 
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            collect <- NULL
            for (i in 1:nthresh) collect <- c(collect, accessD(wd, 
                level = levels[i], boundary = boundary))
            noise.level <- sqrt(dev(collect))
            for (i in 1:nthresh) {
                d <- accessD(wd, level = levels[i], boundary = boundary)
                nd <- length(d)
                neta.d <- (log(nd, base = 2)^(3/2))
                sd2 <- (sum((d/noise.level)^2 - 1)/nd)
                if (verbose == TRUE) {
                  cat("neta.d is ", neta.d, "\nsd2 is ", sd2, 
                    "\n")
                  cat("nd is ", nd, "\n")
                  cat("noise.level ", noise.level, "\n")
                }
                if (force.sure == TRUE || sd2 > neta.d/sqrt(nd)) {
                  if (verbose == TRUE) {
                    cat("SURE: Using SURE\n")
                  }
                  thresh[i] <- sure(d/noise.level)
                }
                else {
                  if (verbose == TRUE) 
                    cat("SURE: (sparse) using sqrt 2log n\n")
                  thresh[i] <- sqrt(2 * log(nd))
                }
                if (verbose == TRUE) 
                  cat("Threshold for level: ", levels[i], " is ", 
                    thresh[i], "\n")
            }
        }
    }
    else if (policy == "BayesThresh") {
        if (alpha < 0) 
            stop("parameter alpha is negative")
        if (beta < 0) 
            stop("parameter beta is negative")
        nthresh <- length(levels)
        nsignal <- rep(0, nthresh)
        noise.level <- sqrt(dev(accessD(wd, level = (nlevelsWT(wd) - 
            1))))
        v <- 2^(-alpha * levels)
        if (is.na(C1)) {
            if (C1.start < 0) 
                stop("C1.start is negative")
            universal <- threshold(wd, policy = "universal", 
                type = "hard", dev = dev, by.level = FALSE, levels = levels)
            sum2 <- rep(0, nthresh)
            for (i in 1:nthresh) {
                dun <- accessD(universal, level = levels[i])
                nsignal[i] <- sum(abs(dun) > 10^-10)
                if (nsignal[i] > 0) 
                  sum2[i] <- sum(dun[abs(dun) > 0]^2)
            }
            if (sum(nsignal) == 0) {
                wd <- nullevels(wd, levelstonu = levels)
                if (verbose == TRUE) 
                  cat("hyperparameters of the prior are: alpha = ", 
                    alpha, "C1 = 0", "beta = ", beta, "C2 = 0\n")
                return(wd)
            }
            else {
                fntoopt <- function(C, nsignal, noise.level, 
                  wd, sum2, v) {
                  ans <- nsignal * (log(noise.level^2 + C^2 * 
                    v) - 2 * log(pnorm((-noise.level * sqrt(2 * 
                    log(2^nlevelsWT(wd))))/sqrt(noise.level^2 + 
                    C^2 * v)))) + sum2/(noise.level^2 + C^2 * 
                    v)
                  sum(ans)
                }
                C1 <- optimize(f = fntoopt, interval = c(0, 50 * 
                  sqrt(C1.start)), nsignal = nsignal, noise.level = noise.level, 
                  wd = wd, sum2 = sum2, v = v)$minimum^2
            }
        }
        if (C1 < 0) 
            stop("parameter C1 is negative")
        tau2 <- C1 * v
        if (is.na(C2)) {
            p <- 2 * pnorm((-noise.level * sqrt(2 * log(2^wd$nlevels)))/sqrt(noise.level^2 + 
                tau2))
            if (beta == 1) 
                C2 <- sum(nsignal/p)/nlevelsWT(wd)
            else C2 <- (1 - 2^(1 - beta))/(1 - 2^((1 - beta) * 
                wd$nlevels)) * sum(nsignal/p)
        }
        if (C2 < 0) 
            stop("parameter C2 is negative")
        if (verbose == TRUE) 
            cat("noise.level is: ", round(noise.level, 4), "\nhyperparameters of the prior are: alpha = ", 
                alpha, "C1 = ", round(C1, 4), "beta = ", beta, 
                "C2 = ", round(C2, 4), "\n")
        if (C1 == 0 | C2 == 0) 
            wd <- nullevels(wd, levelstonu = levels)
        else {
            pr <- pmin(1, C2 * 2^(-beta * levels))
            rat <- tau2/(noise.level^2 + tau2)
            for (i in 1:nthresh) {
                d <- accessD(wd, level = levels[i])
                w <- (1 - pr[i])/pr[i]/sqrt((noise.level^2 * 
                  rat[i])/tau2[i]) * exp((-rat[i] * d^2)/2/noise.level^2)
                z <- 0.5 * (1 + pmin(w, 1))
                d <- sign(d) * pmax(0, rat[i] * abs(d) - noise.level * 
                  sqrt(rat[i]) * qnorm(z))
                wd <- putD(wd, level = levels[i], v = d)
            }
        }
        return(wd)
    }
    else if (policy == "cv") {
        if (verbose == TRUE) 
            cat("Cross-validation policy\n")
        if (by.level == TRUE) 
            stop("Cross-validation policy does not permit by.level\n\t\t\tthresholding (yet)")
        ynoise <- wr(wd)
        thresh <- CWCV(ynoise = ynoise, x = 1:length(ynoise), 
            filter.number = wd$filter$filter.number, family = wd$filter$family, 
            thresh.type = type, tol = cvtol, maxits = cvmaxits, 
            verbose = 0, plot.it = FALSE, ll = min(levels))$xvthresh
        thresh <- rep(thresh, length = nthresh)
    }
    else if (policy == "fdr") {
        if (verbose == TRUE) 
            cat("FDR policy...")
        if (by.level == FALSE) {
            if (verbose == TRUE) 
                cat("All levels at once\n")
            for (i in 1:nthresh) {
                d <- c(d, accessD(wd, level = levels[i], boundary = boundary))
            }
            if (length(value) != 1) 
                stop("Length of value should be 1")
            noise.level <- sqrt(dev(accessD(wd, level = (nlevelsWT(wd) - 
                1))))
            minit <- length(d)
            dinit <- d
            thinit <- qnorm(1 - Q/2) * noise.level
            if (log(n, 2) > 12) 
                ninit <- 3
            else {
                if (log(n, 2) > 10) 
                  ninit <- 2
                else ninit <- 1
            }
            for (k in seq(1, ninit)) {
                dinit1 <- dinit[abs(dinit) >= thinit]
                minit <- length(dinit1)
                if (minit == 0) 
                  thresh <- max(abs(d)) * 1.0001
                else {
                  thinit <- qnorm(1 - (Q * minit)/(2 * n)) * 
                    noise.level
                  minit1 <- length(dinit1[abs(dinit1) >= thinit])
                  if (minit1 == minit || minit1 == 0) 
                    break
                  dinit <- dinit1
                }
            }
            if (noise.level > 0) {
                m <- length(d)
                minit <- length(dinit)
                p <- (2 - 2 * pnorm(abs(dinit)/noise.level))
                index <- order(p)
                j <- seq(1, minit)
                m0 <- max(j[p[index] <= (Q * j)/m])
                if (m0 != "NA" && m0 < minit) 
                  thresh <- abs(dinit[index[m0]])
                else {
                  if (m0 == "NA") 
                    thresh <- max(abs(dinit)) * 1.0001
                  else thresh <- 0
                }
            }
            else thresh <- 0
            thresh <- rep(thresh, length = nthresh)
            if (verbose == TRUE) 
                cat("Global threshold is: ", thresh[1], "\n", 
                  "sigma is: ", noise.level, "\n")
        }
        else {
            if (verbose == TRUE) 
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            for (i in 1:nthresh) {
                d <- accessD(wd, level = levels[i], boundary = boundary)
                m <- length(d)
                noise.level <- sqrt(dev(d))
                thinit <- qnorm(1 - Q/2) * noise.level
                dinit <- d[abs(d) >= thinit]
                minit <- length(dinit)
                if (minit == 0) 
                  thresh[i] <- max(abs(d)) * 1.0001
                else {
                  if (noise.level > 0) {
                    p <- (2 - 2 * pnorm(abs(dinit)/noise.level))
                    index <- order(p)
                    j <- seq(1, minit)
                    m0 <- max(j[p[index] <= (Q * j)/m])
                    if (m0 != "NA" && m0 < minit) 
                      thresh[i] <- abs(dinit[index[m0]])
                    else {
                      if (m0 == "NA") 
                        thresh[i] <- max(abs(dinit)) * 1.0001
                      else thresh[i] <- 0
                    }
                  }
                  else thresh[i] <- 0
                }
                if (verbose == TRUE) 
                  cat("Threshold for level: ", levels[i], "is", 
                    thresh[i], "\n")
            }
        }
    }
    else if (policy == "op1") {
        if (verbose == TRUE) 
            cat("Ogden and Parzen's first policy\n")
        if (by.level == FALSE) 
            stop("Ogden and Parzen's first policy only computes level-dependent policies")
        thresh <- TOthreshda1(ywd = wd, alpha = OP1alpha, verbose = verbose, 
            return.threshold = return.threshold)
        return(thresh)
    }
    else if (policy == "op2") {
        if (verbose == TRUE) 
            cat("Ogden and Parzen's second policy\n")
        if (by.level == FALSE) 
            stop("Ogden and Parzen's second policy only computes level-dependent policies")
        thresh <- TOthreshda2(ywd = wd, alpha = OP1alpha, verbose = verbose, 
            return.threshold = return.threshold)
        return(thresh)
    }
    else if (policy == "manual") {
        if (verbose == TRUE) 
            cat("Manual policy\n")
        thresh <- rep(value, length = nthresh)
        if (length(value) != 1 && length(value) != nthresh) 
            warning("your threshold is not the same length as number of levels")
    }
    else if (policy == "mannum") {
        if (verbose == TRUE) {
            cat("Manual policy using ", value, " of the")
            cat(" largest coefficients\n")
        }
        if (value < 1) {
            stop("Have to select an integer larger than 1 for value")
        }
        else if (value > length(wd$D)) {
            stop(paste("There are only ", length(wd$D), " coefficients, you specified ", 
                value))
        }
        coefs <- wd$D
        scoefs <- sort(abs(coefs))
        scoefs <- min(rev(scoefs)[1:value])
        wd$D[abs(wd$D) < scoefs] <- 0
        return(wd)
    }
    else if (policy == "probability") {
        if (verbose == TRUE) 
            cat("Probability policy...")
        if (by.level == FALSE) {
            if (verbose == TRUE) 
                cat("All levels at once\n")
            for (i in 1:nthresh) d <- c(d, accessD(wd, level = levels[i], 
                boundary = boundary))
            if (length(value) != 1) 
                stop("Length of value should be 1")
            thresh <- rep(quantile(abs(d), prob = value), length = nthresh)
            if (verbose == TRUE) 
                cat("Global threshold is: ", thresh[1], "\n")
        }
        else {
            if (verbose == TRUE) 
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            if (length(value) == 1) 
                value <- rep(value, nthresh)
            if (length(value) != nthresh) 
                stop("Wrong number of probability values")
            for (i in 1:nthresh) {
                d <- accessD(wd, level = levels[i], boundary = boundary)
                thresh[i] <- quantile(abs(d), prob = value[i])
                if (verbose == TRUE) 
                  cat("Threshold for level: ", levels[i], " is ", 
                    thresh[i], "\n")
            }
        }
    }
    if (return.threshold == TRUE) 
        return(thresh)
    for (i in 1:nthresh) {
        d <- accessD(wd, level = levels[i], boundary = boundary)
        if (type == "hard") {
            d[abs(d) <= thresh[i]] <- 0
        }
        else if (type == "soft") {
            d <- (d * (abs(d) - thresh[i]) * (abs(d) > thresh[i]))/abs(d)
            d[is.na(d)] <- 0
        }
        if (verbose == TRUE) 
            cat("Level: ", levels[i], " there are ", sum(d == 
                0), " zeroes\n")
        wd <- putD(wd, level = levels[i], v = d, boundary = boundary)
    }
    wd
}
