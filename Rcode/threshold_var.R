threshold.var <- function (x.w, x.w.v, levels, type = "hard")
{
    d <- NULL
    n <- 2^nlevelsWT(x.w)
    J <- nlevelsWT(x.w)
    nthresh <- length(levels)
    thresh <- list(0)
    for (i in 1:nthresh) {
        d <- accessD(x.w,level=levels[i])
        ind <- (((J-1)-levels[i])*n+1):((J-levels[i])*n)
        noise.level <- sqrt(x.w.v[ind])
        nd <- length(d)
        thresh[[i]] <- sqrt(2*log(nd*log2(nd))) * noise.level
    }
    for (i in 1:nthresh) {
        d <- accessD(x.w,level=levels[i])
        if (type == "hard") {
            d[abs(d) <= thresh[[i]]] <- 0
        }
        else if (type == "soft") {
            d <- (d * (abs(d) - thresh[[i]]) * (abs(d) > thresh[[i]]))/abs(d)
            d[is.na(d)] <- 0
        }
        x.w=putD(x.w,level=levels[i],v=d)
    }
    return(x.w)
}






