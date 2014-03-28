threshold.wd.var <- function (wd, wd.var, levels = 3:(nlevelsWT(wd) - 1), type = "hard",
     policy = "universal", boundary = FALSE, verbose = FALSE, return.threshold = FALSE)
{
    d <- NULL
    n <- 2^nlevelsWT(wd)
    nthresh <- length(levels)
    if (verbose == TRUE) 
        cat("Universal policy...")
    thresh <- list(0)
    for (i in 1:nthresh) {
        d <- accessD(wd, level = levels[i], boundary = boundary)
        noise.level <- sqrt(accessC(wd.var, level = levels[i], boundary = boundary))
        nd <- length(d)
        thresh[[i]] <- sqrt(2 * log(nd)) * noise.level
        if (verbose == TRUE) 
          cat("Threshold for level: ", levels[i], " is ", 
            thresh[[i]], "\n")
    }
        if (return.threshold == TRUE) 
        return(thresh)
    for (i in 1:nthresh) {
        d <- accessD(wd, level = levels[i], boundary = boundary)
        if (type == "hard") {
            d[abs(d) <= thresh[[i]]] <- 0
        }
        else if (type == "soft") {
            d <- (d * (abs(d) - thresh[[i]]) * (abs(d) > thresh[[i]]))/abs(d)
            d[is.na(d)] <- 0
        }
        if (verbose == TRUE) 
            cat("Level: ", levels[i], " there are ", sum(d == 
                0), " zeroes\n")
        wd <- putD(wd, level = levels[i], v = d, boundary = boundary)
    }
    wd
}






