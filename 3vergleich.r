#R <3vergleich.r --no-save

load("results.RData")

library("MASS")
library("robust")


#Ranking nach Standarddeviation
t1 <- with(res[[1]], data.frame(src, id))
for(fktName in names(res)) {
  t2 <- with(res[[fktName]], data.frame(src, id, sd))
  names(t2)[3] <- fktName
  t1 <- merge(t1, t2)
}
##Guttenberg
t2 <- as.double(as.matrix(t1[t1$src=="gut",3:NCOL(t1)]))
t3 <- round(c(0, quantile(t2, probs=(1:6)/7), max(t2+0.01)), 2)
t4 <- apply(t1[t1$src=="gut",3:NCOL(t1)], 1, cut, breaks=t3, labels=F)
t5 <- table(rep(colnames(t1)[3:NCOL(t1)], NCOL(t4)), t4)
t6 <- apply(t(t1[t1$src=="gut",3:NCOL(t1)]), 1, median)
names(t6) <- colnames(t1)[3:NCOL(t1)]
t6 <- sort(t6, decreasing = T)
t5 <- t5[names(t6),]
rownames(t5) <- paste(names(t6), " (", formatC(round(t6, 2), format='f', digits=2), ")", sep="")
colnames(t5) <- paste("<", t3[2:length(t3)], sep="")
pdf("./pic/sdRankGutten.pdf", width=5, height=9.1, pointsize=10)
#pdf("./pic/sdRankGuttenSW.pdf", width=5, height=8.5, pointsize=10)
par(mar=c(0,9,0,0.2))
barplot(t(t5),horiz=TRUE, las=1, col=colorRampPalette(c("#008000", "#00FF00", "#FFFF00", "#FF0000", "#800000"))(NCOL(t5)), space=.3, legend=colnames(t5), args.legend = list(x="bottomright", horiz = T, bty="n"), border=NA)
#barplot(t(t5),horiz=TRUE, las=1, legend=colnames(t5), args.legend = list(x="bottomright", horiz = T, bty="n"), space=.3, col=1, angle=c(135,45), density=c(0,10,10,10,20,20,20))
#rownames(t5) <- NULL
#barplot(t(t5),horiz=TRUE, las=1, legend=colnames(t5), args.legend = list(x="bottomright", horiz = T, bty="n"), space=.3, col=1, angle=c(135,45,135,135,135,45,45), density=c(0,10,10,10,20,20,20), add=T, border=NA)
dev.off()
##BFW
t2 <- as.double(as.matrix(t1[t1$src=="jrm",3:NCOL(t1)]))
t3 <- round(c(0, quantile(t2, probs=(1:6)/7), max(t2+0.01)), 2)
t4 <- apply(t1[t1$src=="jrm",3:NCOL(t1)], 1, cut, breaks=t3, labels=F)
t5 <- table(rep(colnames(t1)[3:NCOL(t1)], NCOL(t4)), t4)
t6 <- apply(t(t1[t1$src=="jrm",3:NCOL(t1)]), 1, median)
names(t6) <- colnames(t1)[3:NCOL(t1)]
t6 <- sort(t6, decreasing = T)
t5 <- t5[names(t6),]
rownames(t5) <- paste(names(t6), " (", formatC(round(t6, 2), format='f', digits=2), ")", sep="")
colnames(t5) <- paste("<", t3[2:length(t3)], sep="")
pdf("./pic/sdRankBFW.pdf", width=5, height=9.1, pointsize=10)
#pdf("./pic/sdRankBFWSW.pdf", width=5, height=8.5, pointsize=10)
par(mar=c(0,9,0,0.2))
barplot(t(t5),horiz=TRUE, las=1, col=colorRampPalette(c("#008000", "#00FF00", "#FFFF00", "#FF0000", "#800000"))(NCOL(t5)), space=.3, legend=colnames(t5), args.legend = list(x="bottomright", horiz = T, bty="n"), border=NA)
#barplot(t(t5),horiz=TRUE, las=1, legend=colnames(t5), args.legend = list(x="bottomright", horiz = T, bty="n"), space=.3, col=1, angle=c(135,45), density=c(0,10,10,10,20,20,20))
#rownames(t5) <- NULL
#barplot(t(t5),horiz=TRUE, las=1, legend=colnames(t5), args.legend = list(x="bottomright", horiz = T, bty="n"), space=.3, col=1, angle=c(135,45,135,135,135,45,45), density=c(0,10,10,10,20,20,20), add=T, border=NA)
dev.off()


#Hoehenunterschiede Hoehenkurve und Beobachtung zu den Zeitpunkten ageSteps
ageSteps <- c(10, 30, 70, 100, 150)
ys <- 1.3
sh <- c(-0.3, -0.15, 0, 0.15, 0.3) * ys
t1 <- data.frame(fun=NA, age=NA, src=NA, id=NA, val=NA)[0,]
for(fktName in names(res)) {
  resCoef <- res[[fktName]]
  for(i in 1:length(ageSteps)) {
    t2 <- data.frame(fun=fktName, age=ageSteps[i], src=resCoef[,"src"], id=resCoef[,"id"], val = resCoef[,paste("h", ageSteps[i], sep="_")] - resCoef[,paste("hEst", ageSteps[i], sep="_")])
    t1 <- rbind(t1, t2[complete.cases(t2),])
  }
}
#Median und Perzentillen (5%, 25%)
t2 <- aggregate(val ~ fun + age, data=t1[t1$src=="gut",], quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
#t2$fun <- factor(t2$fun, levels = sort(levels(t2$fun), decreasing=T))
t2$fun <- factor(t2$fun, levels = sort(unique(t2$fun), decreasing=T))
t2 <- data.frame(t2[,1:2], t2$val)
#t3 <- aggregate(X50. ~ fun, data=t2, function(x) {sum(abs(x))})
#t3 <- data.frame(fun=t3$fun[order(t3$X50.)], line=NROW(t3):1)
t3 <- sort(unlist(data.frame(rbind(by(t2, t2$fun, function(x) {sum(abs(x$X25.)+abs(x$X75.))})))))
t3 <- data.frame(fun=names(t3), line=length(t3):1 * ys)
t2 <- merge(t2, t3)
pdf("./pic/hdiffGutten.pdf", width=8, height=14.6, pointsize=20)
par(mar=c(2,7,0.1,0.1), xaxs="i")
with(t2, plot(X50., line, yaxt="n", xlab="", ylab="", xlim=c(-.75,.75), type="n", ylim=ys*c(2.5,nlevels(t2$fun)-1.5), cex.axis=0.8))
axis(2,at=t3$line,labels=t3$fun, las=1, tick=F, cex.axis=0.8)
abline(h=ys*(0.5 + 1:(nlevels(t2$fun)-1)), lty=3, lwd=1, col="grey")
abline(v=seq(-0.6,0.6,0.2), lty=3, lwd=1, col="grey")
abline(v=0, col="grey")
for(i in 1:length(ageSteps)) {
  with(t2[t2$age==ageSteps[i],], segments(X50., sh[i]+line-0.1, X50., sh[i]+line+0.1, lwd=1))
  with(t2[t2$age==ageSteps[i],], segments(X5., sh[i]+line, X95., lty=3))
  with(t2[t2$age==ageSteps[i],], segments(X25., sh[i]+line, X75., lwd=1))
  with(t2[t2$age==ageSteps[i] & t2$X50. > .75,], points(pmin(.75, X50.), sh[i]+line, pch=19, cex=0.7))
  with(t2[t2$age==ageSteps[i] & t2$X50. < -.75,], points(pmax(-.75, X50.), sh[i]+line, pch=19, cex=0.7))
}
dev.off()
#
t2 <- aggregate(val ~ fun + age, data=t1[t1$src=="jrm",], quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
t2$fun <- factor(t2$fun, levels = sort(unique(t2$fun), decreasing=T))
t2 <- data.frame(t2[,1:2], t2$val)
t3 <- sort(unlist(data.frame(rbind(by(t2, t2$fun, function(x) {sum(abs(x$X25.)+abs(x$X75.))})))))
t3 <- data.frame(fun=names(t3), line=length(t3):1 * ys)
t2 <- merge(t2, t3)
pdf("./pic/hdiffBFW.pdf", width=8, height=14.6, pointsize=20)
par(mar=c(2,7,0.1,0.1), xaxs="i")
with(t2, plot(X50., line, yaxt="n", xlab="", ylab="", xlim=c(-.75,.75), type="n", ylim=ys*c(2.5,nlevels(t2$fun)-1.5), cex.axis=0.8))
axis(2,at=t3$line,labels=t3$fun, las=1, tick=F, cex.axis=0.8)
abline(h=ys*(0.5 + 1:(nlevels(t2$fun)-1)), lty=3, lwd=1, col="grey")
abline(v=seq(-0.6,0.6,0.2), lty=3, lwd=1, col="grey")
abline(v=0, col="grey")
for(i in 1:length(ageSteps)) {
  with(t2[t2$age==ageSteps[i],], segments(X50., sh[i]+line-0.1, X50., sh[i]+line+0.1, lwd=1))
  with(t2[t2$age==ageSteps[i],], segments(X5., sh[i]+line, X95., lty=3))
  with(t2[t2$age==ageSteps[i],], segments(X25., sh[i]+line, X75., lwd=1))
  with(t2[t2$age==ageSteps[i] & t2$X50. > .75,], points(pmin(.75, X50.), sh[i]+line, pch=19, cex=0.7))
  with(t2[t2$age==ageSteps[i] & t2$X50. < -.75,], points(pmax(-.75, X50.), sh[i]+line, pch=19, cex=0.7))
}
dev.off()


#Trend ueber Bonitaet?
ageSteps <- c(10, 30, 70, 100, 150)
ys <- 1.3
sh <- c(-0.3, -0.15, 0, 0.15, 0.3) * ys
t1 <- data.frame(fun=NA, age=NA, src=NA, id=NA, h=NA, hest=NA)[0,]
for(fktName in names(res)) {
  resCoef <- res[[fktName]]
  for(i in 1:length(ageSteps)) {
    t2 <- data.frame(fun=fktName, age=ageSteps[i], src=resCoef[,"src"], id=resCoef[,"id"], h = resCoef[,paste("h", ageSteps[i], sep="_")], hest = resCoef[,paste("hEst", ageSteps[i], sep="_")])
    t1 <- rbind(t1, t2[complete.cases(t2),])
  }
}
#Guttenberg
t2 <- by(t1[t1$src=="gut",], list(t1$fun[t1$src=="gut"], t1$age[t1$src=="gut"]) , function(x) summary(lmRob(x$h-x$hest ~ x$hest))$coefficients[2,c("Estimate","Std. Error")])
t3 <- data.frame(rbind(t2))
t5 <- data.frame(fun=NA, age=NA, val=NA, se=NA)[0,]
for(i in names(t3)) {
  t4 <- as.data.frame(t(data.frame(t3[,i])))
  t4$age <- substring(i, 2)
  t4$fun <- rownames(t4)
  t5 <- rbind(t5, t4[,c(4,3,1,2)])
}
colnames(t5) <- c("fun", "age", "val", "se")
t5$fun <- factor(t5$fun, levels = sort(unique(t5$fun), decreasing=T))
t5$age <- as.numeric(t5$age)
t6 <- sort(unlist(data.frame(rbind(by(t5[t5$age>10,], t5$fun[t5$age>10], function(x) {sum(abs(x$val)+x$se/2)})))))
t6 <- data.frame(fun=names(t6), line=length(t6):1 * ys)
t5 <- merge(t5, t6)
#
pdf("./pic/hdifTrendGut.pdf", width=8, height=14.6, pointsize=20)
par(mar=c(2,7,0.1,0.1), xaxs="i")
#with(t5, plot(val, as.numeric(fun), yaxt="n", xlab="", ylab="", type="n", xlim=range(c(val-se, val+se))))
with(t5, plot(val, as.numeric(fun), yaxt="n", xlab="", ylab="", type="n", xlim=c(-0.51,0.11), ylim=ys*c(2.5,nlevels(t5$fun)-1.5), cex.axis=0.8))
axis(2,at=t6$line,labels=t6$fun, las=1, tick=F, cex.axis=0.8)
abline(h=ys*(0.5 + 1:(nlevels(t5$fun)-1)), lty=3, lwd=1, col="grey")
abline(v=0, col="grey")
for(i in 1:length(ageSteps)) {
  with(t5[t5$age==ageSteps[i],], segments(val, sh[i]+line-0.1, val, sh[i]+line+0.1, lwd=1))
  with(t5[t5$age==ageSteps[i],], segments(val - se, sh[i]+line, val + se, lwd=1))
  with(t5[t5$age==ageSteps[i] & t5$val > .11,], points(pmin(.11, val), sh[i]+line, pch=19, cex=0.7))
  with(t5[t5$age==ageSteps[i] & t5$val < -.51,], points(pmax(-.51, val), sh[i]+line, pch=19, cex=0.7))
}
dev.off()
#BFW
t2 <- by(t1[t1$src=="jrm",], list(t1$fun[t1$src=="jrm"], t1$age[t1$src=="jrm"]) , function(x) summary(lmRob(x$h-x$hest ~ x$hest))$coefficients[2,c("Estimate","Std. Error")])
t3 <- data.frame(rbind(t2))
t5 <- data.frame(fun=NA, age=NA, val=NA, se=NA)[0,]
for(i in names(t3)) {
  t4 <- as.data.frame(t(data.frame(t3[,i])))
  t4$age <- substring(i, 2)
  t4$fun <- rownames(t4)
  t5 <- rbind(t5, t4[,c(4,3,1,2)])
}
colnames(t5) <- c("fun", "age", "val", "se")
t5$fun <- factor(t5$fun, levels = sort(unique(t5$fun), decreasing=T))
t5$age <- as.numeric(t5$age)
t6 <- sort(unlist(data.frame(rbind(by(t5[t5$age>10,], t5$fun[t5$age>10], function(x) {sum(abs(x$val)+x$se/2)})))))
t6 <- data.frame(fun=names(t6), line=length(t6):1 * ys)
t5 <- merge(t5, t6)
#
pdf("./pic/hdifTrendBFW.pdf", width=8, height=14.6, pointsize=20)
par(mar=c(2,7,0.1,0.1), xaxs="i")
with(t5, plot(val, as.numeric(fun), yaxt="n", xlab="", ylab="", type="n", xlim=c(-0.21,0.11), ylim=ys*c(2.5,nlevels(t5$fun)-1.5), cex.axis=0.8))
axis(2,at=t6$line,labels=t6$fun, las=1, tick=F, cex.axis=0.8)
abline(h=ys*(0.5 + 1:(nlevels(t5$fun)-1)), lty=3, lwd=1, col="grey")
abline(v=0, col="grey")
for(i in 1:length(ageSteps)) {
  with(t5[t5$age==ageSteps[i],], segments(val, sh[i]+line-0.1, val, sh[i]+line+0.1, lwd=1))
  with(t5[t5$age==ageSteps[i],], segments(val - se, sh[i]+line, val + se, lwd=1))
  with(t5[t5$age==ageSteps[i] & t5$val > .11,], points(pmin(.11, val), sh[i]+line, pch=19, cex=0.7))
  with(t5[t5$age==ageSteps[i] & t5$val < -.21,], points(pmax(-.21, val), sh[i]+line, pch=19, cex=0.7))
}
dev.off()



#Geschaetzte hoehen ueber dem Alter
ageSteps <- c(10, 30, 70, 100, 150, 300, 800)
ys <- 1.3
sh <- c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3) * ys
t1 <- data.frame(fun=NA, age=NA, src=NA, id=NA, h=NA)[0,]
for(fktName in names(res)) {
  resCoef <- res[[fktName]]
  for(i in 1:length(ageSteps)) {
    t2 <- data.frame(fun=fktName, age=ageSteps[i], src=resCoef[,"src"], id=resCoef[,"id"], h = resCoef[,paste("hEst", ageSteps[i], sep="_")])
    t1 <- rbind(t1, t2[complete.cases(t2),])
  }
}
#Median und Perzentillen (5%, 25%) - Guttenberg
t2 <- aggregate(h ~ fun + age, data=t1[t1$src=="gut",], quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
t2$fun <- factor(t2$fun, levels = sort(unique(t2$fun), decreasing=T))
colnames(t2)[3] <- "val"
t2 <- data.frame(t2[,1:2], t2$val)
t3 <- data.frame(fun = t2$fun[t2$age==800][order(t2$X50.[t2$age==800] - t2$X50.[t2$age==150])], line=1:nlevels(t2$fun) * ys)
t2 <- merge(t2, t3)
#
pdf("./pic/hDevGut.pdf", width=8, height=14.6, pointsize=20)
par(mar=c(2,7,0.1,0.1), xaxs="i")
#with(t2, plot(X50., as.numeric(fun), yaxt="n", xlab="", ylab="", xlim=range(t2[,3:7]), type="n", log="x"))
with(t2, plot(X50., as.numeric(fun), yaxt="n", xlab="", ylab="", xlim=c(0.2,80), type="n", log="x", ylim=ys*c(2.5,nlevels(t2$fun)-1.5), cex.axis=0.8))
#axis(2,at=1:nlevels(t2$fun),labels=levels(t2$fun), las=1, cex.axis=1)
axis(2,at=t3$line,labels=t3$fun, las=1, tick=F, cex.axis=0.8)
abline(h=ys*(0.5 + 1:(nlevels(t2$fun)-1)), lty=3, lwd=1, col="grey")
abline(v=c(1,5,10,15,(2:20)*10), lty=3, col="grey")
for(i in 1:length(ageSteps)) {
  with(t2[t2$age==ageSteps[i],], segments(X50., sh[i]+line-0.1, X50., sh[i]+line+0.1, lwd=1))
  with(t2[t2$age==ageSteps[i],], segments(X5., sh[i]+line, X95., lty=3))
  with(t2[t2$age==ageSteps[i],], segments(X25., sh[i]+line, X75., lwd=1))
  with(t2[t2$age==ageSteps[i] & t2$X50. > 80,], points(pmin(80, X50.), sh[i]+line, pch=19, cex=0.7))
  with(t2[t2$age==ageSteps[i] & t2$X50. < .2,], points(pmax(.2, X50.), sh[i]+line, pch=19, cex=0.7))
}
dev.off()
#BFW
t2 <- aggregate(h ~ fun + age, data=t1[t1$src=="jrm",], quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
t2$fun <- factor(t2$fun, levels = sort(unique(t2$fun), decreasing=T))
colnames(t2)[3] <- "val"
t2 <- data.frame(t2[,1:2], t2$val)
t3 <- data.frame(fun = t2$fun[t2$age==800][order(t2$X50.[t2$age==800] - t2$X50.[t2$age==150])], line=1:nlevels(t2$fun)*ys)
t2 <- merge(t2, t3)
pdf("./pic/hDevBFW.pdf", width=8, height=14.6, pointsize=20)
par(mar=c(2,7,0.1,0.1), xaxs="i")
with(t2, plot(X50., as.numeric(fun), yaxt="n", xlab="", ylab="", xlim=c(0.2,80), type="n", log="x", ylim=ys*c(2.5,nlevels(t2$fun)-1.5), cex.axis=0.8))
axis(2,at=t3$line,labels=t3$fun, las=1, tick=F, cex.axis=0.8)
abline(h=ys*(0.5 + 1:(nlevels(t2$fun)-1)), lty=3, lwd=1, col="grey")
abline(v=c(1,5,10,15,(2:20)*10), lty=3, col="grey")
for(i in 1:length(ageSteps)) {
  with(t2[t2$age==ageSteps[i],], segments(X50., sh[i]+line-0.1, X50., sh[i]+line+0.1, lwd=1))
  with(t2[t2$age==ageSteps[i],], segments(X5., sh[i]+line, X95., lty=3))
  with(t2[t2$age==ageSteps[i],], segments(X25., sh[i]+line, X75., lwd=1))
  with(t2[t2$age==ageSteps[i] & t2$X50. > 80,], points(pmin(80, X50.), sh[i]+line, pch=19, cex=0.7))
  with(t2[t2$age==ageSteps[i] & t2$X50. < .2,], points(pmax(.2, X50.), sh[i]+line, pch=19, cex=0.7))
}
dev.off()
#Negative Hoehenzuwachse?
t2 <- aggregate(h ~ fun + id, data=t1, diff)
t2 <- aggregate(h ~ fun + id, data=t1, function(x) {min(diff(x))})
aggregate(h ~ fun, data=t2, min)
