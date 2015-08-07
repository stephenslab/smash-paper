library(vioplot)

l2norm = function(x) sum(x^2)
mise.sep = function(x, y) 10000 * apply(x - matrix(rep(y, times = 100), nrow = 100, byrow = T), 1, l2norm)/l2norm(y)

mise.BMSM.bur.1 = mise.sep((est.BMSM.bur.1), mu.bur.1)
mise.BMSM.bur.8 = mise.sep((est.BMSM.bur.8), mu.bur.8)


mise.ash.bur.1 = mise.sep(t(est.ash.bur.1), mu.bur.1)
mise.ash.bur.8 = mise.sep(t(est.ash.bur.8), mu.bur.8)

mise.hf.ti.r.4.bur.1 = mise.sep(t(est.hf.ti.r.4.bur.1), mu.bur.1)
mise.hf.ti.r.5.bur.1 = mise.sep(t(est.hf.ti.r.5.bur.1), mu.bur.1)
mise.hf.ti.r.6.bur.1 = mise.sep(t(est.hf.ti.r.6.bur.1), mu.bur.1)
mise.hf.ti.r.7.bur.1 = mise.sep(t(est.hf.ti.r.7.bur.1), mu.bur.1)
mise.hf.ti.r.bur.1 = colMeans(matrix(c(mise.hf.ti.r.4.bur.1, mise.hf.ti.r.5.bur.1, mise.hf.ti.r.6.bur.1, mise.hf.ti.r.7.bur.1), byrow = TRUE, nr = 4))

mise.hf.ti.r.4.bur.8 = mise.sep(t(est.hf.ti.r.4.bur.8), mu.bur.8)
mise.hf.ti.r.5.bur.8 = mise.sep(t(est.hf.ti.r.5.bur.8), mu.bur.8)
mise.hf.ti.r.6.bur.8 = mise.sep(t(est.hf.ti.r.6.bur.8), mu.bur.8)
mise.hf.ti.r.7.bur.8 = mise.sep(t(est.hf.ti.r.7.bur.8), mu.bur.8)
mise.hf.ti.r.bur.8 = colMeans(matrix(c(mise.hf.ti.r.4.bur.8, mise.hf.ti.r.5.bur.8, mise.hf.ti.r.6.bur.8, mise.hf.ti.r.7.bur.8), byrow = TRUE, nr = 4))

pdf("paper/violin_pois_1.pdf", height = 8, width = 12)
vioplot(mise.ash.bur.1, mise.BMSM.bur.1, mise.hf.ti.r.bur.1, names = c("SMASH", "BMSM", "Haar-Fisz"))
title(ylab = "MISE", xlab = "Method")
abline(h = median(mise.ash.bur.1), lty = 3, col = 3)
dev.off()

pdf("paper/violin_pois_8.pdf", height = 8, width = 12)
vioplot(mise.ash.bur.8, mise.BMSM.bur.8, mise.hf.ti.r.bur.8, names = c("SMASH", "BMSM", "Haar-Fisz"))
title(ylab = "MISE", xlab = "Method")
abline(h = median(mise.ash.bur.8), lty = 3, col = 3)
dev.off()