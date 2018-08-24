# This script creates plots showing the intensity functions used to
# simulate the Poisson data.

# SCRIPT PARAMETERS
# -----------------
n = 1024
t = 1:n/n

# GENERATE THE TRUE SIGNALS
# -------------------------
# Some of the intensity functions are defined in signals.R.
source("../code/signals.R")
mu.s   = spike.f(t)
mu.ang = dop.f(t)

sig = ((2 * t + 0.5) * (t <= 0.15)) +
    ((-12 * (t - 0.15) + 0.8) * (t > 0.15 & t <= 0.2)) +
    0.2 * (t > 0.2 & t <= 0.5) + 
    ((6 * (t - 0.5) + 0.2) * (t > 0.5 & t <= 0.6)) +
    ((-10 * (t - 0.6) + 0.8) * (t > 0.6 & t <= 0.65)) +
    ((-0.5 * (t - 0.65) + 0.3) * (t > 0.65 & t <= 0.85)) +
    ((2 * (t - 0.85) + 0.2) * (t > 0.85))
mu.ang = 3/5 * ((5/(max(sig) - min(sig))) * sig - 1.6) - 0.0419569

heavi = 4 * sin(4 * pi * t) - sign(t - 0.3) - sign(0.72 - t)
mu.hs = heavi/sqrt(var(heavi)) * 1 * 2.99/3.366185
mu.hs = mu.hs - min(mu.hs)

I_1 = exp(-(abs(t - 0.2)/0.01)^1.2) * (t <= 0.2) +
      exp(-(abs(t - 0.2)/0.03)^1.2) * (t > 0.2)
I_2 = exp(-(abs(t - 0.3)/0.01)^1.2) * (t <= 0.3) +
      exp(-(abs(t - 0.3)/0.03)^1.2) * (t > 0.3)
I_3 = exp(-(abs(t - 0.4)/0.01)^1.2) * (t <= 0.4) +
      exp(-(abs(t - 0.4)/0.03)^1.2) * (t > 0.4)
mu.bur = 2.99/4.51804 * (4 * I_1 + 3 * I_2 + 4.5 * I_3)

pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
hgt = 2.88/5 * c(4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
mu.cb = rep(0, n)
for (j in 1:length(pos))
  mu.cb = mu.cb + (1 + sign(t - pos[j])) * (hgt[j]/2)
mu.cb[mu.cb < 0] = 0

pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
hgt = 2.97/5 * c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
wth = c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
mu.b = rep(0, n)
for (j in 1:length(pos))
  mu.b = mu.b + hgt[j]/((1 + (abs(t - pos[j])/wth[j]))^4)

# PLOT THE SIGNAL INTENSITIES
# ---------------------------
plot(t,mu.s,xlab = "",ylab = "",type = "l",
     main = "Spikes intensity function")
invisible(readline(prompt = "Press [enter] to continue analysis... "))
plot(t,mu.ang,xlab = "",ylab = "",type = "l",
     main = "Angles intensity function")
invisible(readline(prompt = "Press [enter] to continue analysis... "))
plot(t,mu.hs,xlab = "",ylab = "",type = "l",
     main = "Heavisine intensity function")
invisible(readline(prompt = "Press [enter] to continue analysis... "))
plot(t,mu.bur,xlab = "",ylab = "",type = "l",
     main = "Bursts intensity function")
invisible(readline(prompt = "Press [enter] to continue analysis... "))
plot(t,mu.cb,xlab = "",ylab = "",type = "l",
     main = "Clipped Blocks intensity function")
invisible(readline(prompt = "Press [enter] to continue analysis... "))
plot(t,mu.b,xlab = "",ylab = "",type = "l",
     main = "Bumps intensity function")

