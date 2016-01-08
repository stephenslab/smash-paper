systime = system.time(ashsmooth.gaus(X.s,filter.number=1,family="DaubExPhase"))

tempfile = "ashprofile.txt"
Rprof(tempfile, interval = 0.1)
mu.est.haar<-ashsmooth.gaus(X.s,filter.number=1,family="DaubExPhase")
Rprof(NULL)

summaryRprof(tempfile)
