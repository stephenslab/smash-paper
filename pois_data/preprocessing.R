peaks=read.table("Gm1287peaks_chr1_sorted")
peaks=peaks[order(peaks[,2]),]


all.base=1:peaks[dim(peaks)[1],3]

peaks.start=all.base%in%peaks[,2]
peaks.end=all.base%in%peaks[,3]

plot(peaks.y[880000:1011072],type='s')

plot(peaks.y[2470000:2535536],type='s')


880000:1011072
2470000:2535536