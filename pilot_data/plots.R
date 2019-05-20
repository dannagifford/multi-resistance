rpdf("rampingsdr.pdf", width=8,height=4)
par(mfcol=c(2,3), cex=0.5, las=1)

for(j in c("rif","nal")){
for(i in c("BW100","MUTS100")){
single.temp=subset(single, treatment==i&sample!="Blank B"&antibiotic==j)
beanplot(od~day, what=c(0,1,0,0), single.temp,col=rgb(0,0,0,.2), ylim=c(-0.4,2.1), log="", xlab="Day", ylab="Optical density",main=paste(i,j))
with(single.temp,points(jitter(day), od,col=rgb(0,.4,.4,.6), cex=0.5, pch=19))
axis(side=3, at=1:7, labels=c("1/8","1/4","1/2","1",2,4,8))
rm(single.temp)
abline(h=0.1, col="red", lty="dashed")
}
}

with(subset(BW.stack, group=="B"),{ beanplot(OD.blank~day, what=c(0,1,0,0), col=rgb(0,0,0,.2), ylim=c(-0.4,2.1), ylab="Optical density", xlab="Day", main="BW100 rif+nal"); points(OD.blank~jitter(day), col=rgb(0,.4,.4,.6), cex=0.5, pch=19);abline(h=0.1, col="red", lty="dashed")})
axis(side=3, at=1:7, labels=c("1/8","1/4","1/2","1",2,4,8))
with(subset(MUTS.stack, group=="B"),{ beanplot(OD.blank~day, what=c(0,1,0,0), col=rgb(0,0,0,.2), ylim=c(-0.4,2.1), ylab="Optical density", xlab="Day", main="MUTS100 rif+nal"); points(OD.blank~jitter(day), col=rgb(0,.4,.4,.6), cex=0.5, pch=19);abline(h=0.1, col="red", lty="dashed")})
axis(side=3, at=1:7, labels=c("1/8","1/4","1/2","1",2,4,8))



dev.off()


rmatrix=matrix(nrow=6,ncol=4)
rmatrix[1,]=c(2,2,2,1.5)
rmatrix[2,]=c(1.6,1.8,1.9,1.6)
rmatrix[3,]=c(1.2,1.6,1.8,1.7)
rmatrix[4,]=c(0.8,1.4,1.7,1.8)
rmatrix[5,]=c(0.4,1.2,1.6,1.9)
rmatrix[6,]=c(0.4,1.0,1.5,2)
#rmatrix[4,]=c(1.4,1.4,1.4,2)
#rmatrix[5,]=c(0.8,1.2,1.2,2)
#rmatrix[6,]=c(0,1,1,2)
#rmatrix[7,]=c(0,0.8,0.8,2)
#rmatrix[8,]=c(0,0,0,2)
#rmatrix[9,]=c(0,0,0,2)
#rmatrix[10,]=c(0,0,0,2)


#temp11=doParam(n <-5, param1.name<-"mu1", param1.seq<-unique(sort(matrix((sapply(5:10, function(X) X^-(6:9))),ncol=1))), param2.name="K", param2.seq<-unique(sort(matrix((sapply(3:10, function(X) X^(3:10))),ncol=1))), ,K=1e9,rmatrix=rmatrix, mu=rep(1e-6, times=4), D=1/100, State=c(w=1e8,x=0,y=0,z=0,a=1e3,b=0,c=0,d=0), Time=seq(0,10,length=n), nbot=10, FUN=LotVcompMNEW, Mutator=100)


#temp10=doParam(n <-50, param1.name<-"mu1", param1.seq<-unique(sort(matrix((sapply(5:10, function(X) X^-(6:9))),ncol=1))), param2.name="K", param2.seq<-unique(sort(matrix((sapply(3:10, function(X) X^(3:10))),ncol=1))), ,K=1e9,rmatrix=rmatrix, mu=rep(1e-6, times=4), D=1/100, State=c(w=1e8,x=1,y=0,z=0,a=1e3,b=0,c=0,d=0), Time=seq(0,2,length=n), nbot=100, FUN=LotVcompMNEW, Mutator=100)


temp0=doParam(n <-5, param1.name<-"mu1", param1.seq<-10^(seq(-10,-7,.5)), param2.name="K", param2.seq<-10^(seq(1,8,.5)), K=1e9,rmatrix=rmatrix, mu=rep(1e-6, times=4), D=1/100, State=c(w=.99,x=0,y=0,z=0,a=.1,b=0,c=0,d=0), Time=seq(0,7,length=n), nbot=5, FUN=LotVcompMNEW, Mutator=100)
temp0.melt=melt(temp0[-length(temp0)])
temp0.rsh=cbind(setNames(temp0.melt[temp0.melt$Var2=="time",c("value","L2","L1")], c("time","mu","LK")),sapply(letters[c(1:4,23:26)], function(i)setNames(temp0.melt[temp0.melt$Var2==i,"value"],i)))
temp0.rsh2=melt(temp0.rsh, id.vars=c("time","mu","LK"), measure.vars=c(letters[c(1:4,23:26)]))
names(temp0.rsh2)[4]="clone"
names(temp0.rsh2)[5]="N"
for(i in c(1,2,3,5)){temp0.rsh2[,i]=as.numeric(temp0.rsh2[,i])}
for(i in unique(temp0.rsh2$time)){print(levelplot(N/LK~log10(mu) + log10(LK) | clone , data = subset(temp0.rsh2, time==i), pretty=F, cuts=1000, col.regions = colorRampPalette(rev(brewer.pal(n=8, name="Spectral")))(1001), interpolate=F, useRaster=F, region=T, layout=c(4,2), at=seq(from=0, to=1, by=0.001))) }



#temp0=doParam(n <-5, param1.name<-"mu1", param1.seq<-c(1e-6,1e-7,1e-8), param2.name="K", param2.seq<-10^(5:7), K=1e9,rmatrix=rmatrix, mu=rep(1e-6, times=4), D=1/100, State=c(w=1e3,x=0,y=0,z=0), Time=seq(0,7,length=n), nbot=10, FUN=LotVcompNEW, Mutator=100)
tempM=lapply(temp0[-length(temp0)],function(Y) sapply(Y, function(X) X[nrow(X),2:ncol(X)]))

tempM=lapply(temp0[-length(temp0)],function(Y) sapply(Y, function(X) X[nrow(X),2:ncol(X)]))



temp3M=melt(tempM, value.name="N",varnames=c("clone","mu"), measure.vars="K", level="K")
temp3M$LK=as.numeric(temp3M$LK)
levelplot(N/LK~log10(mu) + log10(LK) | clone , data = temp3M, pretty=F, cuts=1000, col.regions = colorRampPalette(rev(brewer.pal(n=8, name="Spectral")))(1001), interpolate=F, useRaster=F, region=T, layout=c(4,2)) 



#################### proportion of each clone from ODE
for(i in as.character(10^(5:9))){

test9=setNames((lapply(names(solutions), function(Y) t(sapply(solutions[Y][[1]][[i]], function(X)X[nrow(X),2:9])))), nm=names(solutions))

test99=setNames(melt(test9, measure.vars=letters[c(1:4,23:26)]), c("mu","clone","N","p"))
test99$p=as.numeric(gsub("[a-z]","",test99$p))
for(i in c(1,3,4))test99[,i]=as.numeric(test99[,i])
test999=dcast(data=test99, p + mu ~clone , value.var="N")

real=lapply(ps, function(PS){doParam(n <-5, param1.name<-"mu1", param1.seq<-10^(seq(-10,-7,.5)), param2.name="K", param2.seq<-10^(seq(1,8,.5)), K=1e9,rmatrix=rmatrix, mu=rep(1e-6, times=4), D=1/100, State=c(w=1-PS,x=0,y=0,z=0,a=PS,b=0,c=0,d=0), Time=seq(0,7,length=n), nbot=5, FUN=LotVcompMNEW, Mutator=100)})
names(real)=ps

getFinal=function(PS){
mat1=setNames((lapply(names(PS), function(Y) t(sapply(PS[Y][[1]][[i]], function(X)X[nrow(X),2:9])))), nm=names(PS))

mat2=setNames(melt(mat1, measure.vars=letters[c(1:4,23:26)]), c("mu","clone","N","p"))

#mat2$p=as.numeric(gsub("[a-z]","",test99$p))
for(i in c(1,3,4))mat2[,i]=as.numeric(mat2[,i])

mat3=dcast(data=mat2, p + mu ~clone , value.var="N")
return(mat3)

}


par(las=1, mfrow=c(2,2))

with(test999, plot(y=(a+w)/(a+b+c+d+w+x+y+z), x=mu, log="x",pch=as.numeric(as.factor(p)), type="p", xlab="log10 mutation rate", ylab="Proportion sensitive",xaxt="n", ylim=c(0,1), main="Sensitive (A+W)"))
#with(test999, points(y=(a/(a+b+c+d+w+x+y+z)), x=mu, col="red", log="x",cex=0.5,pch=as.numeric(as.factor(p))))


legend(x="topright", legend=c(rev(unique(test999$p))[1:4],"0.0001"), bty="n", pch=as.numeric(as.factor(rev(unique(test999$p)))), title="p")
axis(side=1, at=10^(-10:-7), labels=-c(10:7))

with(test999, plot(y=(b+x)/(a+b+c+d+w+x+y+z), x=mu, log="x",pch=as.numeric(as.factor(p)), type="p", xlab="log10 mutation rate", ylab="Proportion mutant 1",xaxt="n", ylim=c(0,1), main="Mutant 1 (B+X)"))
#with(test999, points(y=(b/(a+b+c+d+w+x+y+z)), x=mu, col="red", log="x",cex=0.5,pch=as.numeric(as.factor(p))))


legend(x="topleft", legend=c(rev(unique(test999$p))[1:4],"0.0001"), bty="n", pch=as.numeric(as.factor(rev(unique(test999$p)))), title="p")
axis(side=1, at=10^(-10:-7), labels=-c(10:7))

with(test999, plot(y=(c+y)/(a+b+c+d+w+x+y+z), x=mu, log="x",pch=as.numeric(as.factor(p)), type="p", xlab="log10 mutation rate", ylab="Proportion mutant 2",xaxt="n", ylim=c(0,1), main="Mutant 2 (C+Y)"))
axis(side=1, at=10^(-10:-7), labels=-c(10:7))
legend(x="topright", legend=c(rev(unique(test999$p))[1:4],"0.0001"), bty="n", pch=as.numeric(as.factor(rev(unique(test999$p)))), title="p")
#with(test999, points(y=(c/(a+b+c+d+w+x+y+z)), x=mu, col="red", log="x",cex=0.5,pch=as.numeric(as.factor(p))))

with(test999, plot(y=(d+z)/(a+b+c+d+w+x+y+z), x=mu, log="x",pch=as.numeric(as.factor(p)), type="p", xlab="log10 mutation rate", ylab="Proportion double mutant",xaxt="n", ylim=c(0,1), main="Double-mutant (D+Z)"))
axis(side=1, at=10^(-10:-7), labels=-c(10:7))
legend(x="topleft", legend=c(rev(unique(test999$p))[1:4],"0.0001"), bty="n", pch=as.numeric(as.factor(rev(unique(test999$p)))), title="p")
#with(test999, points(y=(d/(a+b+c+d+w+x+y+z)), x=mu, col="red", log="x",cex=0.5,pch=as.numeric(as.factor(p))))

}
with(test999, plot(y=(x/(a+b+c+d+w+x+y+z)), x=mu, col="blue", log="x",cex=0.5,pch=as.numeric(as.factor(p)), type="n", ylim=c(0,1)))
for(i in unique(test999$p)){
with(subset(test999, p==i), points(y=(d/(a+b+c+d+w+x+y+z)), x=mu, col="blue", log="x",cex=0.5,pch=as.numeric(as.factor(p==i))))
with(subset(test999, p==i), points(y=(z/(a+b+c+d+w+x+y+z)), x=mu, col="red", log="x",cex=0.5,pch=as.numeric(as.factor(p[p==i]))))

}
with(test999, points(y=(a/(a+b+c+d+w+x+y+z)), x=mu, col="red", log="x",cex=0.5,pch=as.numeric(as.factor(p))))
with(test999, points(y=(b/(a+b+c+d+w+x+y+z)), x=mu, col="red", log="x",cex=0.5,pch=as.numeric(as.factor(p))))
with(test999, points(y=(c/(a+b+c+d+w+x+y+z)), x=mu, col="red", log="x",cex=0.5,pch=as.numeric(as.factor(p))))
with(test999, points(y=(d/(a+b+c+d+w+x+y+z)), x=mu, col="red", log="x",cex=0.5,pch=as.numeric(as.factor(p))))

real=lapply(ps, function(PS){doParam(n <-5, param1.name<-"mu1", param1.seq<-10^(seq(-10,-7,.5)), param2.name="K", param2.seq<-10^(seq(1,8,.5)), K=1e9,rmatrix=rmatrix, mu=rep(1e-6, times=4), D=1/100, State=c(w=1-PS,x=0,y=0,z=0,a=PS,b=0,c=0,d=0), Time=seq(0,7,length=n), nbot=5, FUN=LotVcompMNEW, Mutator=100)})
