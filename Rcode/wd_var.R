wd.var <- function (data, type = "station", filter.number = 8, family = "DaubLeAsymm", 
    bc = "periodic", verbose = FALSE, min.scale = 0, precond = TRUE) 
{
    if (verbose == TRUE) 
        cat("wd: Argument checking...")
    if (!is.atomic(data)) 
        stop("Data is not atomic")
    DataLength <- length(data)
    nlevels <- nlevelsWT(data)
    if (is.na(nlevels)) 
        stop("Data length is not power of two")
    if (type != "wavelet" && type != "station") 
        stop("Unknown type of wavelet decomposition")
    if (type == "station" && bc != "periodic") 
        stop("Can only do periodic boundary conditions with station")
    if (verbose == TRUE) 
        cat("...done\nFilter...")
    if (bc != "interval") 
        filter <- filter.select(filter.number, family)
        filter$H <- filter.select(filter.number,family)$H^2
    if (verbose == TRUE) 
        cat("...selected\nFirst/last database...")
    fl.dbase <- first.last(LengthH = length(filter$H), DataLength = DataLength, 
        type = type, bc = bc)
    if (bc == "interval") {
        ans <- wd.int(data = data, preferred.filter.number = filter.number, 
            min.scale = min.scale, precond = precond)
        fl.dbase <- first.last(LengthH = length(filter$H), DataLength = DataLength, 
            type = type, bc = bc, current.scale = min.scale)
        filter <- list(name = paste("CDV", filter.number, sep = ""), 
            family = "CDV", filter.number = filter.number)
        l <- list(transformed.vector = ans$transformed.vector, 
            current.scale = ans$current.scale, filters.used = ans$filters.used, 
            preconditioned = ans$preconditioned, date = ans$date, 
            nlevels = IsPowerOfTwo(length(ans$transformed.vector)), 
            fl.dbase = fl.dbase, type = type, bc = bc, filter = filter)
        class(l) <- "wd"
        return(l)
    }
    dtsp <- tsp(data)
    C <- rep(0, fl.dbase$ntotal)
    C[1:DataLength] <- data
    if (verbose == TRUE) 
        error <- 1
    else error <- 0
    if (verbose == TRUE) 
        cat("built\n")
    if (verbose == TRUE) 
        cat("Decomposing...\n")
    nbc <- switch(bc, periodic = 1, symmetric = 2)
    if (is.null(nbc)) 
        stop("Unknown boundary condition")
    ntype <- switch(type, wavelet = 1, station = 2)
    if (is.null(filter$G)) {
        wavelet.decomposition <- .C("wavedecomp", C = as.double(C), 
            D = as.double(rep(0, fl.dbase$ntotal.d)), H = as.double(filter$H), 
            LengthH = as.integer(length(filter$H)), nlevels = as.integer(nlevels), 
            firstC = as.integer(fl.dbase$first.last.c[, 1]), 
            lastC = as.integer(fl.dbase$first.last.c[, 2]), offsetC = as.integer(fl.dbase$first.last.c[, 
                3]), firstD = as.integer(fl.dbase$first.last.d[, 
                1]), lastD = as.integer(fl.dbase$first.last.d[, 
                2]), offsetD = as.integer(fl.dbase$first.last.d[, 
                3]), ntype = as.integer(ntype), nbc = as.integer(nbc), 
            error = as.integer(error), PACKAGE = "wavethresh")
    }
    else {
        wavelet.decomposition <- .C("comwd", CR = as.double(Re(C)), 
            CI = as.double(Im(C)), LengthC = as.integer(fl.dbase$ntotal), 
            DR = as.double(rep(0, fl.dbase$ntotal.d)), DI = as.double(rep(0, 
                fl.dbase$ntotal.d)), LengthD = as.integer(fl.dbase$ntotal.d), 
            HR = as.double(Re(filter$H)), HI = as.double(-Im(filter$H)), 
            GR = as.double(Re(filter$G)), GI = as.double(-Im(filter$G)), 
            LengthH = as.integer(length(filter$H)), nlevels = as.integer(nlevels), 
            firstC = as.integer(fl.dbase$first.last.c[, 1]), 
            lastC = as.integer(fl.dbase$first.last.c[, 2]), offsetC = as.integer(fl.dbase$first.last.c[, 
                3]), firstD = as.integer(fl.dbase$first.last.d[, 
                1]), lastD = as.integer(fl.dbase$first.last.d[, 
                2]), offsetD = as.integer(fl.dbase$first.last.d[, 
                3]), ntype = as.integer(ntype), nbc = as.integer(nbc), 
            error = as.integer(error), PACKAGE = "wavethresh")
    }
    if (verbose == TRUE) 
        cat("done\n")
    error <- wavelet.decomposition$error
    if (error != 0) {
        cat("Error ", error, " occured in wavedecomp\n")
        stop("Error")
    }
    if (is.null(filter$G)) {
        l <- list(C = wavelet.decomposition$C, D = wavelet.decomposition$D, 
            nlevels = nlevelsWT(wavelet.decomposition), fl.dbase = fl.dbase, 
            filter = filter, type = type, bc = bc, date = date())
    }
    else {
        l <- list(C = complex(real = wavelet.decomposition$CR, 
            imaginary = wavelet.decomposition$CI), D = complex(real = wavelet.decomposition$DR, 
            imaginary = wavelet.decomposition$DI), nlevels = nlevelsWT(wavelet.decomposition), 
            fl.dbase = fl.dbase, filter = filter, type = type, 
            bc = bc, date = date())
    }
    class(l) <- "wd"
    if (!is.null(dtsp)) 
        tsp(l) <- dtsp
    l
}
