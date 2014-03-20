library(wavethresh)
source("ash/bayesmooth.R")
source("ash/bayesmooth_test.R")


tt=wd(1:16,8,"DaubLeAsymm","station")
tt=wd.var.s8(c(2,rep(1,15)),"station")

accessC(tt,4)

accessC(tt,3)
sum(hf^2*accessC(tt,4))
sum((c(hf[15:16]^2,hf[1:14]^2)*accessC(tt,4)))
sum(hf^2*rshift(accessC(tt,4)))
sum((c(hf[15:16]^2,hf[1:14]^2)*rshift(accessC(tt,4))))



sum(hf*accessC(tt,2))
sum((c(hf[15:16],hf[1:14])*accessC(tt,2)))
sum(hf*rshift(accessC(tt,2)))
sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2))))

accessD(tt,1)
2*sum((hf*accessC(tt,2))[c(1,3,5,7,9,11,13,15)])-sum(hf*accessC(tt,2))
2*sum((c(hf[15:16],hf[1:14])*accessC(tt,2))[c(1,3,5,7,9,11,13,15)])-sum((c(hf[15:16],hf[1:14])*accessC(tt,2)))
2*sum((hf*rshift(accessC(tt,2)))[c(1,3,5,7,9,11,13,15)])-sum(hf*rshift(accessC(tt,2)))
2*sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2)))[c(1,3,5,7,9,11,13,15)])-sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2))))


accessC(tt,0)
sum(hf*c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2)))))
sum(hf*c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2))))))
sum(hf*rshift(c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2))))))
sum(hf*rshift(c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2)))))))



accessD(tt,0)
2*sum((hf*c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2)))))[c(1,3,5,7,9,11,13,15)])-sum(hf*c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2)))))
2*sum((c(hf[15:16],hf[1:14])*c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2))))))
2*sum((hf*rshift(c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*rshift(c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2))))))
2*sum((c(hf[15:16],hf[1:14])*rshift(c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2)))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*rshift(c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2)))))))


2*sum((hf*c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2)))))[c(1,3,5,7,9,11,13,15)])-sum(hf*c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2)))))+5
-(2*sum((hf*c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2)))))[c(1,3,5,7,9,11,13,15)])-sum(hf*c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2))))))+5
2*sum((c(hf[15:16],hf[1:14])*c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2))))))+5
-(2*sum((c(hf[15:16],hf[1:14])*c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2)))))))+5

2*sum((hf*rshift(c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*rshift(c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2))))))+5
-(2*sum((hf*rshift(c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*rshift(c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2)))))))+5
2*sum((c(hf[15:16],hf[1:14])*rshift(c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2)))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*rshift(c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2)))))))+5
-(2*sum((c(hf[15:16],hf[1:14])*rshift(c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2)))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*rshift(c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2))))))))+5

0.5*(2*sum((hf*c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2)))))[c(1,3,5,7,9,11,13,15)])-sum(hf*c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2)))))+5+lshift(2*sum((hf*rshift(c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*rshift(c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2))))))+5))
0.5*(-(2*sum((hf*c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2)))))[c(1,3,5,7,9,11,13,15)])-sum(hf*c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2))))))+5+lshift(-(2*sum((hf*rshift(c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*rshift(c(sum(hf*accessC(tt,2)),sum((c(hf[15:16],hf[1:14])*accessC(tt,2)))))))+5))
0.5*(2*sum((c(hf[15:16],hf[1:14])*c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2))))))+5+lshift(2*sum((c(hf[15:16],hf[1:14])*rshift(c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2)))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*rshift(c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2)))))))+5))
0.5*(-(2*sum((c(hf[15:16],hf[1:14])*c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2)))))))+5+lshift(-(2*sum((c(hf[15:16],hf[1:14])*rshift(c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2)))))))[c(1,3,5,7,9,11,13,15)])-sum(hf*rshift(c(sum(hf*rshift(accessC(tt,2))),sum((c(hf[15:16],hf[1:14])*rshift(accessC(tt,2))))))))+5))



,name="Daub cmpct on least asymm N=8,family="DaubLeAsymm",filter.number=8
