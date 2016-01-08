---
title: "ASH profiling in SMASH"
output: html_document
---

Tests the three optimization methods in ASH: "mixEM", "cxxMixSquarem" and "mixIP"


```r
library(SMASH)
```

```
## Error in library(SMASH): there is no package called 'SMASH'
```

```r
load("data.Robj")

wd = smash:::titable(X.s)$difftable
wv = smash:::titable(sigma.t^2)$sumtable

scale = 10

time.em = matrix(0, nr = 100, nc = 2)
time.squarem = matrix(0, nr = 100, nc = 2)
time.ip = matrix(0, nr = 100, nc = 2)
for(i in 1:100){
  time.em[i, ] = system.time(ash(wd[scale, ], sqrt(wv[scale, ]), prior = "nullbiased", optmethod = "mixEM", pointmass = TRUE, nullcheck = TRUE, mixsd = NULL, mixcompdist = "normal", gridmult = 2, df = NULL))[1:2]
  time.squarem[i, ] = system.time(ash(wd[scale, ], sqrt(wv[scale, ]), prior = "nullbiased", optmethod = "cxxMixSquarem", pointmass = TRUE, nullcheck = TRUE, mixsd = NULL, mixcompdist = "normal", gridmult = 2, df = NULL))[1:2]
  time.ip[i, ] = system.time(ash(wd[scale, ], sqrt(wv[scale, ]), prior = "nullbiased", optmethod = "mixIP", pointmass = TRUE, nullcheck = TRUE, mixsd = NULL, mixcompdist = "normal", gridmult = 2, df = NULL))[1:2]
}

summary(time.em)
```

```
##        V1               V2        
##  Min.   :0.0500   Min.   :0.0000  
##  1st Qu.:0.1100   1st Qu.:0.0000  
##  Median :0.1200   Median :0.0000  
##  Mean   :0.1157   Mean   :0.0058  
##  3rd Qu.:0.1300   3rd Qu.:0.0100  
##  Max.   :0.1500   Max.   :0.0400
```

```r
summary(time.squarem)
```

```
##        V1               V2        
##  Min.   :0.0300   Min.   :0.0000  
##  1st Qu.:0.0500   1st Qu.:0.0000  
##  Median :0.0600   Median :0.0000  
##  Mean   :0.0669   Mean   :0.0044  
##  3rd Qu.:0.0900   3rd Qu.:0.0000  
##  Max.   :0.1400   Max.   :0.0400
```

```r
summary(time.ip)
```

```
##        V1               V2        
##  Min.   :0.0700   Min.   :0.0000  
##  1st Qu.:0.1000   1st Qu.:0.0000  
##  Median :0.1100   Median :0.0000  
##  Mean   :0.1174   Mean   :0.0061  
##  3rd Qu.:0.1200   3rd Qu.:0.0100  
##  Max.   :0.3000   Max.   :0.0500
```

```r
scale = 2

time.em = matrix(0, nr = 100, nc = 2)
time.squarem = matrix(0, nr = 100, nc = 2)
time.ip = matrix(0, nr = 100, nc = 2)
for(i in 1:100){
  time.em[i, ] = system.time(ash(wd[scale, ], sqrt(wv[scale, ]), prior = "nullbiased", optmethod = "mixEM", pointmass = TRUE, nullcheck = TRUE, mixsd = NULL, mixcompdist = "normal", gridmult = 2, df = NULL))[1:2]
  time.squarem[i, ] = system.time(ash(wd[scale, ], sqrt(wv[scale, ]), prior = "nullbiased", optmethod = "cxxMixSquarem", pointmass = TRUE, nullcheck = TRUE, mixsd = NULL, mixcompdist = "normal", gridmult = 2, df = NULL))[1:2]
  time.ip[i, ] = system.time(ash(wd[scale, ], sqrt(wv[scale, ]), prior = "nullbiased", optmethod = "mixIP", pointmass = TRUE, nullcheck = TRUE, mixsd = NULL, mixcompdist = "normal", gridmult = 2, df = NULL))[1:2]
}

summary(time.em)
```

```
##        V1               V2        
##  Min.   :0.0900   Min.   :0.0000  
##  1st Qu.:0.1675   1st Qu.:0.0000  
##  Median :0.1800   Median :0.0000  
##  Mean   :0.1866   Mean   :0.0087  
##  3rd Qu.:0.2000   3rd Qu.:0.0200  
##  Max.   :0.3900   Max.   :0.0400
```

```r
summary(time.squarem)
```

```
##        V1              V2        
##  Min.   :0.050   Min.   :0.0000  
##  1st Qu.:0.060   1st Qu.:0.0000  
##  Median :0.070   Median :0.0000  
##  Mean   :0.072   Mean   :0.0013  
##  3rd Qu.:0.080   3rd Qu.:0.0000  
##  Max.   :0.160   Max.   :0.0400
```

```r
summary(time.ip)
```

```
##        V1               V2        
##  Min.   :0.2200   Min.   :0.0000  
##  1st Qu.:0.2900   1st Qu.:0.0000  
##  Median :0.3100   Median :0.0100  
##  Mean   :0.3146   Mean   :0.0156  
##  3rd Qu.:0.3300   3rd Qu.:0.0200  
##  Max.   :0.5000   Max.   :0.0700
```
