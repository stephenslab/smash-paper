threshold.haar <- function (vdtable, vtable, levels, type = "hard", policy = "universal")
{
    wmean <- vdtable[-1,]/2
    d <- NULL
    n <- dim(vdtable)[2]
    nthresh <- length(levels)
    thresh <- list(0)
    for (i in 1:nthresh) {
        d <- vdtable[levels[i]+2,]
        noise.level <- sqrt(vtable[levels[i]+2,])
        nd <- length(d)
        thresh[[i]] <- sqrt(2*log(nd*log2(nd))) * noise.level
    }
    for (i in 1:nthresh) {
        d <- vdtable[levels[i]+2,]
        if (type == "hard") {
            d[abs(d) <= thresh[[i]]] <- 0
        }
        else if (type == "soft") {
            d <- (d * (abs(d) - thresh[[i]]) * (abs(d) > thresh[[i]]))/abs(d)
            d[is.na(d)] <- 0
        }
        wmean[levels[i]+1,] <- d/2
    }
    return(wmean)
}






