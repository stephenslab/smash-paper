### Functions for Gaussian simulations
n = 1024
t = 1:n/n


spike.f = function(x) (0.75 * exp(-500 * (x - 0.23)^2) + 1.5 * exp(-2000 * (x - 0.33)^2) + 3 * exp(-8000 * (x - 0.47)^2) + 
                         2.25 * exp(-16000 * (x - 0.69)^2) + 0.5 * exp(-32000 * (x - 0.83)^2))
mu.s = spike.f(t)

pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
hgt = 2.97/5 * c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
wth = c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
mu.b = rep(0, n)
for (j in 1:length(pos)) {
  mu.b = mu.b + hgt[j]/((1 + (abs(t - pos[j])/wth[j]))^4)
}


dop.f = function(x) sqrt(x * (1 - x)) * sin((2 * pi * 1.05)/(x + 0.05))
mu.dop = dop.f(t)
mu.dop = 3/(max(mu.dop) - min(mu.dop)) * (mu.dop - min(mu.dop))
mu.dop.var = 10 * dop.f(t)
mu.dop.var = mu.dop.var - min(mu.dop.var)



sig = ((2 * t + 0.5) * (t <= 0.15)) + ((-12 * (t - 0.15) + 0.8) * (t > 0.15 & t <= 0.2)) + 0.2 * (t > 0.2 & t <= 0.5) + 
  ((6 * (t - 0.5) + 0.2) * (t > 0.5 & t <= 0.6)) + ((-10 * (t - 0.6) + 0.8) * (t > 0.6 & t <= 0.65)) + ((-0.5 * (t - 
                                                                                                                   0.65) + 0.3) * (t > 0.65 & t <= 0.85)) + ((2 * (t - 0.85) + 0.2) * (t > 0.85))
mu.ang = 3/5 * ((5/(max(sig) - min(sig))) * sig - 1.6) - 0.0419569



heavi = 4 * sin(4 * pi * t) - sign(t - 0.3) - sign(0.72 - t)
mu.hs = heavi/sqrt(var(heavi)) * 1 * 2.99/3.366185
mu.hs = mu.hs - min(mu.hs)


pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
hgt = 2.88/5 * c(4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
mu.blk = rep(0, n)
for (j in 1:length(pos)) {
  mu.blk = mu.blk + (1 + sign(t - pos[j])) * (hgt[j]/2)
}

mu.cblk = mu.blk
mu.cblk[mu.cblk < 0] = 0

mu.blip = (0.32 + 0.6 * t + 0.3 * exp(-100 * (t - 0.3)^2)) * (t >= 0 & t <= 0.8) + (-0.28 + 0.6 * t + 0.3 * exp(-100 * 
                                                                                                                  (t - 1.3)^2)) * (t > 0.8 & t <= 1)


mu.cor = 623.87 * t^3 * (1 - 2 * t) * (t >= 0 & t <= 0.5) + 187.161 * (0.125 - t^3) * t^4 * (t > 0.5 & t <= 0.8) + 3708.470441 * 
  (t - 1)^3 * (t > 0.8 & t <= 1)
mu.cor = (0.6/(max(mu.cor) - min(mu.cor))) * mu.cor
mu.cor = mu.cor - min(mu.cor) + 0.2

mu.sp = (1 + mu.s)/5
mu.bump = (1 + mu.b)/5
mu.blk = 0.2 + 0.6 * (mu.blk - min(mu.blk))/max(mu.blk - min(mu.blk))
mu.ang = (1 + mu.ang)/5
mu.dop = (1 + mu.dop)/5
mu.blip = mu.blip
mu.cor = mu.cor


var1 = rep(1, n)
var2 = (1e-02 + 4 * (exp(-550 * (t - 0.2)^2) + exp(-200 * (t - 0.5)^2) + exp(-950 * (t - 0.8)^2)))
var3 = (1e-02 + 2 * mu.dop.var)
var4 = 1e-02 + mu.b
var5 = 1e-02 + 1 * (mu.cblk - min(mu.cblk))/max(mu.cblk)

var1 = var1/2
var2 = var2/max(var2)
var3 = var3/max(var3)
var4 = var4/max(var4)
var5 = var5/max(var5)



pdf("paper/supplementary/gaus_sp.pdf", width = 6, height = 4)
plot(t, mu.sp, xlab = "", ylab = "", type = "l", main = "Mean function - Spikes")
dev.off()

pdf("paper/supplementary/gaus_bump.pdf", width = 6, height = 4)
plot(t, mu.bump, xlab = "", ylab = "", type = "l", main = "Mean function - Bumps")
dev.off()

pdf("paper/supplementary/gaus_blk.pdf", width = 6, height = 4)
plot(t, mu.blk, xlab = "", ylab = "", type = "l", main = "Mean function - Blocks")
dev.off()

pdf("paper/supplementary/gaus_ang.pdf", width = 6, height = 4)
plot(t, mu.ang, xlab = "", ylab = "", type = "l", main = "Mean function - Angles")
dev.off()

pdf("paper/supplementary/gaus_dop.pdf", width = 6, height = 4)
plot(t, mu.dop, xlab = "", ylab = "", type = "l", main = "Mean function - Doppler")
dev.off()

pdf("paper/supplementary/gaus_blip.pdf", width = 6, height = 4)
plot(t, mu.blip, xlab = "", ylab = "", type = "l", main = "Mean function - Blip")
dev.off()

pdf("paper/supplementary/gaus_cor.pdf", width = 6, height = 4)
plot(t, mu.cor, xlab = "", ylab = "", type = "l", main = "Mean function - Corner")
dev.off()



pdf("paper/supplementary/gaus_var_cons.pdf", width = 6, height = 4)
plot(t, var1, xlab = "", ylab = "", ylim = c(0, 1), type = "l", main = "Variance function - Constant")
dev.off()

pdf("paper/supplementary/gaus_var_texp.pdf", width = 6, height = 4)
plot(t, var2, xlab = "", ylab = "", type = "l", main = "Variance function - Triple Exponential")
dev.off()

pdf("paper/supplementary/gaus_var_dop.pdf", width = 6, height = 4)
plot(t, var3, xlab = "", ylab = "", type = "l", main = "Variance function - Doppler")
dev.off()

pdf("paper/supplementary/gaus_var_bump.pdf", width = 6, height = 4)
plot(t, var4, xlab = "", ylab = "", type = "l", main = "Variance function - Bumps")
dev.off()

pdf("paper/supplementary/gaus_var_cblk.pdf", width = 6, height = 4)
plot(t, var5, xlab = "", ylab = "", type = "l", main = "Variance function - Clipped Blocks")
dev.off()






# Functions for Poisson simulations
n = 1024
t = 1:n/n

spike.f = function(x) (0.75 * exp(-500 * (x - 0.23)^2) + 1.5 * exp(-2000 * (x - 0.33)^2) + 3 * exp(-8000 * (x - 0.47)^2) + 
                         2.25 * exp(-16000 * (x - 0.69)^2) + 0.5 * exp(-32000 * (x - 0.83)^2))
mu.s = spike.f(t)

dop.f = function(x) sqrt(x * (1 - x)) * sin((2 * pi * 1.05)/(x + 0.05))
mu.ang = dop.f(t)

sig = ((2 * t + 0.5) * (t <= 0.15)) + ((-12 * (t - 0.15) + 0.8) * (t > 0.15 & t <= 0.2)) + 0.2 * (t > 0.2 & t <= 0.5) + 
  ((6 * (t - 0.5) + 0.2) * (t > 0.5 & t <= 0.6)) + ((-10 * (t - 0.6) + 0.8) * (t > 0.6 & t <= 0.65)) + ((-0.5 * (t - 
                                                                                                                   0.65) + 0.3) * (t > 0.65 & t <= 0.85)) + ((2 * (t - 0.85) + 0.2) * (t > 0.85))
mu.ang = 3/5 * ((5/(max(sig) - min(sig))) * sig - 1.6) - 0.0419569


heavi = 4 * sin(4 * pi * t) - sign(t - 0.3) - sign(0.72 - t)
mu.hs = heavi/sqrt(var(heavi)) * 1 * 2.99/3.366185
mu.hs = mu.hs - min(mu.hs)

I_1 = exp(-(abs(t - 0.2)/0.01)^1.2) * (t <= 0.2) + exp(-(abs(t - 0.2)/0.03)^1.2) * (t > 0.2)
I_2 = exp(-(abs(t - 0.3)/0.01)^1.2) * (t <= 0.3) + exp(-(abs(t - 0.3)/0.03)^1.2) * (t > 0.3)
I_3 = exp(-(abs(t - 0.4)/0.01)^1.2) * (t <= 0.4) + exp(-(abs(t - 0.4)/0.03)^1.2) * (t > 0.4)
mu.bur = 2.99/4.51804 * (4 * I_1 + 3 * I_2 + 4.5 * I_3)


pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
hgt = 2.88/5 * c(4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
mu.cb = rep(0, n)
for (j in 1:length(pos)) {
  mu.cb = mu.cb + (1 + sign(t - pos[j])) * (hgt[j]/2)
}
mu.cb[mu.cb < 0] = 0

pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
hgt = 2.97/5 * c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
wth = c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
mu.b = rep(0, n)
for (j in 1:length(pos)) {
  mu.b = mu.b + hgt[j]/((1 + (abs(t - pos[j])/wth[j]))^4)
}


pdf("paper/supplementary/pois_sp.pdf", width = 6, height = 4)
plot(t, mu.s, xlab = "", ylab = "", type = "l", main = "Intensity function - Spikes")
dev.off()

pdf("paper/supplementary/pois_ang.pdf", width = 6, height = 4)
plot(t, mu.ang, xlab = "", ylab = "", type = "l", main = "Intensity function - Angles")
dev.off()

pdf("paper/supplementary/pois_hs.pdf", width = 6, height = 4)
plot(t, mu.hs, xlab = "", ylab = "", type = "l", main = "Intensity function - Heavisine")
dev.off()

pdf("paper/supplementary/pois_bur.pdf", width = 6, height = 4)
plot(t, mu.bur, xlab = "", ylab = "", type = "l", main = "Intensity function - Bursts")
dev.off()

pdf("paper/supplementary/pois_cblk.pdf", width = 6, height = 4)
plot(t, mu.cb, xlab = "", ylab = "", type = "l", main = "Intensity function - Clipped Blocks")
dev.off()

pdf("paper/supplementary/pois_bump.pdf", width = 6, height = 4)
plot(t, mu.b, xlab = "", ylab = "", type = "l", main = "Intensity function - Bumps")
dev.off() 