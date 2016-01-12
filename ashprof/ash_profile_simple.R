library(smash)

time_em = matrix(0, nr = 100, nc = 2)
time_squarem = matrix(0, nr = 100, nc = 2)
time_ip = matrix(0, nr = 100, nc = 2)

sd = 2

for(i in 1:10){
  b = rnorm(100000, sd = sd)
  time_em = system.time(ash(b, 1, optmethod = "mixEM"))[1:2]
  time_squarem = system.time(ash(b, 1, optmethod = "cxxMixSquarem"))[1:2]
  time_ip = system.time(ash(b, 1, optmethod = "mixIP"))[1:2]
}

summary(time_em)
summary(time_squarem)
summary(time_ip)



time_em = matrix(0, nr = 100, nc = 2)
time_squarem = matrix(0, nr = 100, nc = 2)
time_ip = matrix(0, nr = 100, nc = 2)

sd = 0.01

for(i in 1:10){
  b = rnorm(100000, sd = sd)
  time_em = system.time(ash(b, 1, optmethod = "mixEM"))[1:2]
  time_squarem = system.time(ash(b, 1, optmethod = "cxxMixSquarem"))[1:2]
  time_ip = system.time(ash(b, 1, optmethod = "mixIP"))[1:2]
}

summary(time_em)
summary(time_squarem)
summary(time_ip)