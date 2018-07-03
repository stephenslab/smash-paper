
# Load the old results.
load("~/git/smash-paper/output/gaus_dscr.RData")
out0 <- res

# Load the new results.
load("res.RData")
out1 <- res

relerr <- function (x, y)
  abs(x - y)/abs(x)

iter <- 0
e    <- rep(0,1400)
methods <- unique(out0$method)
varfn.short <- c("cons", "texp", "dop", "bump", "cblk")
for (m in methods) {
  for (i in c("sp", "bump", "blk", "ang", "dop", "blip", "cor")) {
    for (j in c(1,3)) {
      for (k in 1:5) {
        k0 <- paste0("v",k)
        k1 <- varfn.short[k]
        rows0 <- with(out0,which(method   == m &
                                 scenario == paste(i,j,k0,sep = ".")))
        rows1 <- with(out1,which(method == m &
                                 scenario == paste(i,j,k1,sep = ".")))
        x       <- head(out0$mise[rows0],n = 3)
        y       <- out1$mise[rows1]
        iter    <- iter + 1
        e[iter] <- max(relerr(x,y))
        if (e[iter] > 0.05) {
          print(c(m,i,j,k1,e[iter]))
        }
      }
    }
  }
}
