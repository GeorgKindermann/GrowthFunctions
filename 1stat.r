#R <1stat.r --no-save

#load("results.RData")

library(minpack.lm)
library(nls2)

x <- read.csv("cleanData.csv")
ageSteps <- c(0, 0.0001, 10, 30, 70, 100, 150, 300, 800)
res <- list()

#Hoehenentwicklung
t1 <- unique(data.frame(id=x$id, src=x$src,oh=NA,o50=NA))  #Baumidentifiers
for(i in 1:NROW(t1)) {
  t2 <- x[x$id == t1$id[i],]  #Data of single tree
  a <- splinefun(t2$t, t2$h, method="monoH.FC")
  t1$oh[i] <- a(100)
  t1$o50[i] <- a(50)
}
t1 <- t1[order(t1$src, t1$oh*0.6 + t1$o50),]
t1$dy <- 0
t1$dy[t1$src=="jrm"] <- 29
pdf("./pic/ageHeight.pdf", width=5, height=8.5, pointsize=10)
par(mar=c(4,4,0.2,0.2))
plot(x$t, x$h, type="n", xlab="Age [Years]", ylab="Height [m]", ylim=c(5,233))
abline(h=0:24*10, col="lightgrey")
for(i in 1:NROW(t1)) {
  t2 <- x[x$id == t1$id[i],]  #Data of single tree
  t2 <- t2[order(t2$t),]
  lines(t2$t, i+t2$h-1+t1$dy[i], lty=i)
}
text(0, 40, labels="Guttenberg", srt=90, adj=c(0.5,-0.3))
text(0, 160, labels="BFW", srt=90, adj=c(0.5,-0.3))
dev.off()

#Datenuebersicht
t1 <- with(x, unique(data.frame(id, src)))
table(t1$src)

summary(with(x[x$src=="gut",], tapply(t, id, max)))
summary(with(x[x$src=="jrm",], tapply(t, id, max)))

t1 <- x
t1$src <- as.character(t1$src)
t1$id <- as.character(t1$id)
#Hoehe im alter 50
summary(unlist(lapply(split(t1[t1$src=="gut",], t1$id[t1$src=="gut"]), function(d) splinefun(d$t, d$h, method="monoH.FC")(50))))
summary(unlist(lapply(split(t1[t1$src=="jrm",], t1$id[t1$src=="jrm"]), function(d) splinefun(d$t, d$h, method="monoH.FC")(50))))

#Hoehe im alter 100
summary(unlist(lapply(split(t1[t1$src=="gut",], t1$id[t1$src=="gut"]), function(d) splinefun(d$t, d$h, method="monoH.FC")(100))))
summary(unlist(lapply(split(t1[t1$src=="jrm",], t1$id[t1$src=="jrm"]), function(d) splinefun(d$t, d$h, method="monoH.FC")(100))))
 

#Linear
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=6*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0", "hEst", "si"), rep(ageSteps,each=6), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  a <- lm(h ~ t -1, data=t2)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t4 <- data.frame(t = ageSteps[j], h = t3(ageSteps[j]))
      a <- lm(h ~ t - 1, data=t4)
      resCoef[i,paste(c("c0"), ageSteps[j], sep="_")] <- coef(a)
      #resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Linear"]] <- resCoef
summary(resCoef)


#Parabel
fn <- function(t, c1, tGiv, hGiv) {
  tMul <- hGiv / tGiv + c1*tGiv^2
  tMul * t + c1*t^2
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0", "c1", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  a0 <- lm(h ~ t + I(t^2) -1, data=t2)
  a <- lm(h ~ t + I(t^2) -1, data=t2)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0", "c1")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,tGiv,hGiv), data=t2, start=list(c1=coef(a0)[2]), control=nls.lm.control(maxiter=999), trace=F)
      #t4 <- t3(ageSteps[j]) / (coef(a)[1]*ageSteps[j]^2)
      #resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      #resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Parabola"]] <- resCoef
summary(resCoef)


#Poly3
fn <- function(t, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / tGiv + c1*tGiv^2 + c2*tGiv^3
  tMul * t + c1*t^2 + c2*t^3
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  a0 <- lm(h ~ t + I(t^2) + I(t^3) -1, data=t2)
  a <- a0
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0", "c1", "c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=list(c1=coef(a0)[2],c2=coef(a0)[3]), control=nls.lm.control(maxiter=999), trace=F)
      #t4 <- t3(ageSteps[j]) / (coef(a)[1]*ageSteps[j]^2)
      #resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      #resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Poly3"]] <- resCoef
summary(resCoef)


#Poly4
fn <- function(t, c1,c2,c3, tGiv, hGiv) {
  tMul <- hGiv / tGiv + c1*tGiv^2 + c2*tGiv^3 + c3*tGiv^4
  tMul * t + c1*t^2 + c2*t^3 + c3*t^4
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=9*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=9), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  a0 <- lm(h ~ t + I(t^2) + I(t^3) + I(t^4) -1, data=t2)
  a <- a0
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0", "c1", "c2", "c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=list(c1=coef(a0)[2],c2=coef(a0)[3],c3=coef(a0)[4]), control=nls.lm.control(maxiter=999), trace=F)
      #t4 <- t3(ageSteps[j]) / (coef(a)[1]*ageSteps[j]^2)
      #resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      #resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Poly4"]] <- resCoef
summary(resCoef)


#PolyAllo
fn <- function(t, c1,c2,c3, tGiv, hGiv) {
  tMul <- hGiv / tGiv + c1*tGiv^2 + c2*tGiv^c3
  tMul * t + c1*t^2 + c2*t^c3
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=9*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=9), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  a0 <- lm(h ~ t + I(t^2) + I(t^3) -1, data=t2)
  a0 <- as.numeric(coef(a0))
  a0 <- nlsLM(h ~ c0*t + c1*t^2 + c2*t^c3, data=t2, start=list(c0=a0[1],c1=a0[2],c2=a0[3],c3=3), control=nls.lm.control(maxiter=999), trace=F)
  a <- a0
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0", "c1", "c2", "c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 1) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999), trace=F)
      #t4 <- t3(ageSteps[j]) / (coef(a)[1]*ageSteps[j]^2)
      #resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      #resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["PolyAllo"]] <- resCoef
summary(resCoef)


#Allometric
fn <- function(t, c1, tGiv, hGiv) {
  tMul <- hGiv / (tGiv^c1)
  tMul * (t^c1)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=6*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1", "hEst", "si"), rep(ageSteps,each=6), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- lm(log(h) ~ log(t), data=t2[t2$h>0,])
  a0 <- nls2(h ~ c0*t^c1, data=t2, start=expand.grid(c0=exp(coef(a0)[1]),c1=coef(a0)[2]), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*t^c1, data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,tGiv,hGiv), data=t2, start=coef(a0)[2], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / (ageSteps[j]^coef(a)[1])
      resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Allometric"]] <- resCoef
summary(resCoef)

#Siven
fn <- function(t, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv / ((tGiv/c1)^(c2*tGiv^c3))
  tMul * (t/c1)^(c2*t^c3)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0*(t/c1)^(c2*t^c3), data=t2, start=expand.grid(c0=seq(10,50,10),c1=300,c2=4, c3=-0.5), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*(t/c1)^(c2*t^c3), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999, ftol=0.001), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2", "c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  #a0 <- a
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 1) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / ((ageSteps[j]/coef(a)[1])^(coef(a)[2]*ageSteps[j]^coef(a)[3]))
      resCoef[i,paste(c("c0","c1","c2","c3"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Siven"]] <- resCoef
summary(resCoef)


#Gram
fn <- function(t, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / (tGiv^c1*exp(c2*tGiv))
  tMul * (t^c1*exp(c2*t))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- lm(log(h) ~ log(t) + t, data=t2[t2$h>0,])
  a0 <- nls2(h ~ c0*t^c1*exp(c2*t), data=t2, start=expand.grid(c0=exp(coef(a0)[1]),c1=coef(a0)[2],c2=coef(a0)[3]), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*t^c1*exp(c2*t), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / (ageSteps[j]^coef(a)[1]*exp(coef(a)[2]*ageSteps[j]))
      resCoef[i,paste(c("c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Gram"]] <- resCoef
summary(resCoef)




#Terazaki
#m0 <- nls(H ~ c0*exp(c1*T^-1), start=list(c0=45, c1=-45), data=x)
fn <- function(t, c1, tGiv, hGiv) {
  tMul <- hGiv / exp(c1*tGiv^-1)
  tMul * exp(c1*t^-1)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=6*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1", "hEst", "si"), rep(ageSteps,each=6), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  a0 <- nls2(h ~ c0*exp(c1*t^-1), data=t2, start=expand.grid(c0=seq(10,100,10),c1=seq(-100,-20,10)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*exp(c1*t^-1), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 1) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,tGiv,hGiv), data=t2, start=coef(a0)[2], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / (coef(a)[2]*exp(coef(a)[1]*ageSteps[j]^-1))
      resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Terazaki"]] <- resCoef
summary(resCoef)


#TerazakiE1
fn <- function(t, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / exp(c1*tGiv^-1 + c2*tGiv^-0.5)
  tMul * exp(c1*t^-1 + c2*t^-0.5)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  a0 <- nls2(h ~ c0*exp(c1*t^-1 + c2*t^-0.5), data=t2[t2$t>0,], start=expand.grid(c0=seq(10,100,10),c1=seq(-100,-20,10),c2=0), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*exp(c1*t^-1 + c2*t^-0.5), data=t2[t2$t>0,], start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 1) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2[t2$t>0,], start=coef(a0)[2:3], control=nls.lm.control(maxiter=999), trace=F)
      #t4 <- t3(ageSteps[j]) / (coef(a)[2]*exp(coef(a)[1]*ageSteps[j]^-1))
      #resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["TerazakiE1"]] <- resCoef
summary(resCoef)


#TerazakiE2
fn <- function(t, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv / exp(c1*tGiv^-1 + c2*tGiv^c3)
  tMul * exp(c1*t^-1 + c2*t^c3)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  a0 <- nls2(h ~ c0*exp(c1*t^-1 + c2*t^c3), data=t2[t2$t>0,], start=expand.grid(c0=seq(10,100,10),c1=seq(-100,-20,10),c2=0,c3=-0.5), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*exp(c1*t^-1 + c2*t^c3), data=t2[t2$t>0,], start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2[t2$t>0,], start=coef(a0)[2:4], control=nls.lm.control(maxiter=999), trace=F)
      #t4 <- t3(ageSteps[j]) / (coef(a)[2]*exp(coef(a)[1]*ageSteps[j]^-1))
      #resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["TerazakiE2"]] <- resCoef
summary(resCoef)


#Korf
#m0 <- nls(H ~ c0*exp(c1*T^c2), start=list(c0=35, c1=-35, c2=-1), data=x)
fn <- function(t, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / exp(c1*tGiv^c2)
  tMul * exp(c1*t^c2)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0*exp(c1*t^c2), data=t2, start=expand.grid(c0=seq(10,80,10),c1=seq(-40,-10,5),c2=seq(-1.4,-0.2,.2)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*exp(c1*t^c2), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / exp(coef(a)[1]*ageSteps[j]^coef(a)[2])
      resCoef[i,paste(c("c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Korf"]] <- resCoef
summary(resCoef)



#Gompertz -> geht nicht durch 0!!!
#m0 <- nls(H ~ c0*exp(c1*exp(c2*t)), start=list(c0=35, c1=-5, c2=-0.05), data=x)
fn <- function(t, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / exp(c1*exp(c2*tGiv))
  tMul * exp(c1*exp(c2*t))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0*exp(c1*exp(c2*t)), data=t2, start=expand.grid(c0=seq(10,60,10),c1=seq(-7,-2,0.5),c2=seq(-0.055,-0.015,.005)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*exp(c1*exp(c2*t)), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / exp(coef(a)[1]*exp(coef(a)[2]*ageSteps[j]))
      resCoef[i,paste(c("c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Gompertz"]] <- resCoef
summary(resCoef)


#Sloboda -> Muss nicht automatisch durch 0 gehen!!!
fn <- function(t, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv / exp(c1*exp(c2*tGiv^c3))
  tMul * exp(c1*exp(c2*t^c3))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  #a0 <- nls2(h ~ c0*exp(c1*exp(c2*t^c3)), data=t2, start=expand.grid(c0=seq(20,60,10),c1=seq(-65,-5,10),c2=seq(-2.1,-0.1,.4),c3=seq(.1,1.1,.2)), trace=F, algorithm = "brute-force")
  a0 <- nls2(h ~ c0*exp(c1*exp(c2*t^c3)), data=t2, start=expand.grid(c0=seq(10,80,10),c1=c(-200,-100,-30,-20,-10,-7,-3),c2=c(-4,-3,-2,-1,-0.5,-0.3,-0.2,-0.1,-0.05),c3=seq(.1,1.1,.2)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*exp(c1*exp(c2*t^c3)), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / exp(coef(a)[1]*exp(coef(a)[2]*ageSteps[j]^coef(a)[3]))
      resCoef[i,paste(c("c0","c1","c2","c3"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Sloboda"]] <- resCoef
summary(resCoef)



#Korsun
#m0 <- nls(H ~ c0*exp(c1*log(t)+c2*log(t)^2), start=list(c0=0.0001, c1=5, c2=-0.5), data=x)
fn <- function(t, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / exp(c1*log(tGiv)+c2*log(tGiv)^2)
  tMul * exp(c1*log(t)+c2*log(t)^2)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0*exp(c1*log(t)+c2*log(t)^2), data=t2, start=expand.grid(c0=seq(0.0001,0.03,0.005),c1=seq(1,8,1),c2=seq(-0.8,-0.1,0.1)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*exp(c1*log(t)+c2*log(t)^2), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)
      t4 <- t3(ageSteps[j]) / exp(coef(a)[1]*log(ageSteps[j])+coef(a)[2]*log(ageSteps[j])^2)
      resCoef[i,paste(c("c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Korsun"]] <- resCoef
summary(resCoef)

#KorsunE1
fn <- function(t, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv / exp(c1*log(tGiv)+c2*log(tGiv)^2+c3*log(tGiv)^4)
  tMul * exp(c1*log(t)+c2*log(t)^2+c3*log(t)^4)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a00 <- lm(log(h) ~ log(t) + I(log(t)^2) + I(log(t)^4), data=t2[t2$t>0,])
  a0 <- nls2(h ~ c0*exp(c1*log(t)+c2*log(t)^2+c3*log(t)^4), data=t2[t2$t>0,], start=expand.grid(c0=exp(coef(a00)[1]),c1=coef(a00)[2],c2=coef(a00)[3],c3=coef(a00)[4]), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*exp(c1*log(t)+c2*log(t)^2+c3*log(t)^4), data=t2[t2$t>0,], start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  #resCoef$sd[i] <- summary(a)$sigma
  tmpHest <- coef(a)[1]*exp(coef(a)[2]*log(t2$t) + coef(a)[3]*log(t2$t)^2 + coef(a)[4]*log(t2$t)^4)
  tmpHest[!is.finite(tmpHest)] <- 0
  resCoef$sd[i] <- sqrt(sum((t2$h - tmpHest)^2)/(NROW(t2)-4))
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2[t2$t>0,], start=coef(a0)[2:4], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / exp(coef(a)[1]*log(ageSteps[j])+coef(a)[2]*log(ageSteps[j])^2+coef(a)[3]*log(ageSteps[j])^4)
      resCoef[i,paste(c("c0","c1","c2","c3"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      #resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["KorsunE1"]] <- resCoef
summary(resCoef)


#KosrunE2
fn <- function(t, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv / exp(c1*log(tGiv)+c2*log(tGiv)^c3)
  tMul * exp(c1*log(t)+c2*log(t)^c3)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0*exp(c1*log(t)+c2*log(t)^c3), data=t2[t2$t>0,], start=expand.grid(c0=seq(0.0001,0.03,0.005),c1=seq(1,8,1),c2=seq(-0.8,-0.1,0.1),c3=seq(1,3,0.5)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*exp(c1*log(t)+c2*log(t)^c3), data=t2[t2$t>0,], start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  #resCoef$sd[i] <- summary(a)$sigma
  tmpHest <- coef(a)[1]*exp(coef(a)[2]*log(t2$t) + coef(a)[3]*log(t2$t)^coef(a)[4])
  tmpHest[!is.finite(tmpHest)] <- 0
  resCoef$sd[i] <- sqrt(sum((t2$h - tmpHest)^2)/(NROW(t2)-4))
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2[t2$t>0,], start=coef(a0)[2:4], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / exp(coef(a)[1]*log(ageSteps[j])+coef(a)[2]*log(ageSteps[j])^coef(a)[3])
      resCoef[i,paste(c("c0","c1","c2","c3"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      #resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["KorsunE2"]] <- resCoef
summary(resCoef)


#Weber
fn <- function(t, c1, tGiv, hGiv) {
  tMul <- hGiv / ((1 - exp(c1*tGiv)))
  tMul * (1 - exp(c1*t))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=6*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1", "hEst", "si"), rep(ageSteps,each=6), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0*(1 - exp(c1*t)), data=t2, start=expand.grid(c0=seq(10,80,10),c1=seq(-0.035,-0.005,0.005)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*(1 - exp(c1*t)), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,tGiv,hGiv), data=t2, start=coef(a0)[2], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / ((1 - exp(coef(a)[1]*ageSteps[j])))
      resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Weber"]] <- resCoef
summary(resCoef)


#WeberE1
fn <- function(t, c1,c2, tGiv, hGiv) {
  tMul <- hGiv / ((1 - exp(c1*tGiv))*(1 - exp(c2*tGiv)))
  tMul * (1 - exp(c1*t)) * (1 - exp(c2*t))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  if(i==162) {
    a <- nlsLM(h ~ c0*(1 - exp(c1*t))*(1 - exp(c2*t)), data=t2, start=list(c0=125,c1=-0.0049,c2=-0.0048), lower=c(2,-Inf,-Inf), control=nls.lm.control(maxiter=999), trace=F) } else {
      a <- nlsLM(h ~ c0*(1 - exp(c1*t))*(1 - exp(c2*t)), data=t2, start=list(c0=40,c1=-0.02,c2=-0.01), lower=c(2,-Inf,-Inf), control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)
    }
  a0 <- a
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      if((i==127 && j==6) || (i==149)) {
        a <- nls(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.control(maxiter=999, warnOnly=T), trace=F, algo="port") } else {
          a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)
          resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
        }
      #t4 <- t3(ageSteps[j]) / ((1 - exp(coef(a)[1]*ageSteps[j])))
      #resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["WeberE1"]] <- resCoef
summary(resCoef)


#WeberE2
fn <- function(t, c1,c2,c3, tGiv, hGiv) {
  tMul <- hGiv / ((1 - exp(c1*tGiv))*(1 - exp(c2*tGiv))*(1 - exp(c3*tGiv)))
  tMul * (1 - exp(c1*t)) * (1 - exp(c2*t)) * (1 - exp(c3*t))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
i <- 47
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a <- tryCatch({nlsLM(h ~ c0*(1 - exp(c1*t))*(1 - exp(c2*t))*(1 - exp(c3*t)), data=t2, start=list(c0=40,c1=-0.02,c2=-0.01,c3=-0.001), lower=c(2,-Inf,-Inf,-Inf), control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)},
                error =function(cond) {nlsLM(h ~ c0*(1 - exp(c1*t))*(1 - exp(c2*t))*(1 - exp(c3*t)), data=t2, start=list(c0=45,c1=-0.08,c2=-0.07,c3=-0.01), lower=c(2,-Inf,-Inf,-Inf), control=nls.lm.control(maxiter=999, ftol=0.001), trace=F)})
  a0 <- a
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- tryCatch({nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)},
                    error =function(cond) {nls(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.control(maxiter=0, warnOnly=T), algo="port", trace=F)})
      #t4 <- t3(ageSteps[j]) / ((1 - exp(coef(a)[1]*ageSteps[j])))
      #resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- tryCatch({summary(a)$sigma}, error =function(cond) {NA})
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["WeberE2"]] <- resCoef
summary(resCoef)


#WeberE3
fn <- function(t, c1,c2,c3, tGiv, hGiv) {
  tMul <- hGiv / ((1 - exp(c1*tGiv))*(1 - exp(c2*tGiv))^c3)
  tMul * (1 - exp(c1*t)) * (1 - exp(c2*t))^c3
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
i <- 47
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a <- nlsLM(h ~ c0*(1 - exp(c1*t))*(1 - exp(c2*t))^c3, data=t2, start=list(c0=40,c1=-0.02,c2=-0.01,c3=1), lower=c(2,-Inf,-Inf,-Inf), control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)
  a0 <- a
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)
      #t4 <- t3(ageSteps[j]) / ((1 - exp(coef(a)[1]*ageSteps[j])))
      #resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["WeberE3"]] <- resCoef
summary(resCoef)



#Mitscherlich Richards
#Kurve bei Regression durch vorgegebenen Punkt (tGiv, hGiv) zwingen
fn <- function(t, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / ((1 - exp(c1*tGiv))^c2)
  tMul * (1 - exp(c1*t))^c2
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0*(1 - exp(c1*t))^c2, data=t2, start=expand.grid(c0=seq(10,80,10),c1=seq(-0.035,-0.005,0.005),c2=seq(1,3,.25)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*(1 - exp(c1*t))^c2, data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  #t3 <- with(t2, approxfun(t,h))
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  #Possibilities to get sigma
  #summary(a)$sigma
  #sqrt(sum(resid(a)^2)/(NROW(t2)-3))
  #sqrt(sum((t2$h - coef(a)[1]*(1 - exp(coef(a)[2]*t2$t))^coef(a)[3])^2)/(NROW(t2)-3))
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.lm.control(maxiter=999), trace=F)
    #a0 <- nls2(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=expand.grid(c1=seq(-0.035,-0.005,0.005),c2=seq(1,3,.25)), trace=F, algorithm = "brute-force")
    #a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / ((1 - exp(coef(a)[1]*ageSteps[j]))^coef(a)[2])
      resCoef[i,paste(c("c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Mitscherlich"]] <- resCoef
summary(resCoef)


#Fischer
fn <- function(t, c1, c3, tGiv, hGiv) {
  tMul <- hGiv / ((1 - exp(c1*tGiv^c3)))
  tMul * (1 - exp(c1*t^c3))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c3", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0*(1 - exp(c1*t^c3)), data=t2, start=expand.grid(c0=seq(10,80,10),c1=seq(-0.035,-0.005,0.005),c3=seq(0.5,2,0.5)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*(1 - exp(c1*t^c3)), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / ((1 - exp(coef(a)[1]*ageSteps[j]^coef(a)[2])))
      resCoef[i,paste(c("c0","c1","c3"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Fischer"]] <- resCoef
summary(resCoef)



#Todorovic
fn <- function(t, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv / ((1 - exp(c1*tGiv^c3))^c2)
  tMul * (1 - exp(c1*t^c3))^c2
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0*(1-exp(c1*t^c3))^c2, data=t2, start=expand.grid(c0=seq(10,80,10),c1=c(-4,-1.5,-0.5,-0.1,-0.01,-0.001),c2=c(0.5,seq(1,3,.25),5,8,15,20,50,1000),c3=seq(0.1,2,.2)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*(1-exp(c1*t^c3))^c2, data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / ((1 - exp(coef(a)[1]*ageSteps[j]^coef(a)[3]))^coef(a)[2])
      resCoef[i,paste(c("c0","c1","c2","c3"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Todorovic"]] <- resCoef
summary(resCoef)


#Hossfeld3
#m0 <- nls(H ~ c0*t/(c1 + c1*log(t) + t), start=list(c0=35, c1=8000, c2=-2), data=x)
fn <- function(t, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / ((tGiv) /(c1 + c2*log(tGiv) + tGiv))
  tMul * ((t) /(c1 + c2*log(t) + t))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ (c0 * t) /(c1 + c2*log(t) + t), data=t2, start=expand.grid(c0=seq(5,80,10),c1=seq(200,500,80),c2=seq(-100,-10,10)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ (c0 * t) /(c1 + c2*log(t) + t), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999, ftol=0.001), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.lm.control(maxiter=999, ftol=0.001), trace=F)
      t4 <- t3(ageSteps[j]) / (ageSteps[j]/(coef(a)[1] + coef(a)[2]*log(ageSteps[j]) + ageSteps[j]))
      resCoef[i,paste(c("c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Hossfeld3"]] <- resCoef
summary(resCoef)


#Hossfeld1
#m0 <- nls(H ~ c0/(1+c1/t+c2/t^2), start=list(c0=35, c1=-20, c2=3000), data=x[x$T>0,])
fn <- function(t, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / (1/(1+c1/tGiv+c2/tGiv^2))
  tMul * 1/(1+c1/t+c2/t^2)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0/(1+c1/t+c2/t^2), data=t2, start=expand.grid(c0=seq(20,80,10),c1=seq(10,100,10),c2=seq(2000,6000,1000)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0/(1+c1/t+c2/t^2), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999, ftol=0.001), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.lm.control(maxiter=999, ftol=0.001), trace=F)
      t4 <- t3(ageSteps[j]) / (1/(1+coef(a)[1]/ageSteps[j]+coef(a)[2]/ageSteps[j]^2))
      resCoef[i,paste(c("c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Hossfeld1"]] <- resCoef
summary(resCoef)


#Hos1E1
fn <- function(t, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv / (1/(1+c1/tGiv+c2*tGiv^(-2) + c3*tGiv^(-3)))
  tMul * 1/(1+c1/t+c2*t^(-2) + c3*t^(-3))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a00 <- lm(I(1/h) ~ I(1/t) + I(1/t^2) + I(1/t^3), data=t2[t2$t>0,])
  #a0 <- nls2(h ~ c0/(1+c1/t+c2*t^(-2) + c3*t^(-3)), data=t2[t2$t>0,], start=expand.grid(c0=1/coef(a00)[1],c1=coef(a00)[2]/coef(a00)[1],c2=coef(a00)[3]/coef(a00)[1],c3=coef(a00)[4]/coef(a00)[1]), trace=F, algorithm = "brute-force")
  #a0 <- nls2(h ~ c0/(1+c1/t+c2*t^(-2) + c3*t^(-3)), data=t2[t2$t>0,], start=expand.grid(c0=seq(10,80,10),c1=seq(5,200,20),c2=seq(100,6000,1000),c3=seq(-3,-1,0.5)), trace=F, algorithm = "brute-force")
  a0 <- nls2(h ~ c0/(1+c1/t+c2*t^(-2) + c3*t^(-3)), data=t2[t2$t>0,], start=expand.grid(c0=c(1/coef(a00)[1],seq(10,80,10)),c1=c(coef(a00)[2]/coef(a00)[1],seq(5,200,20)),c2=c(coef(a00)[3]/coef(a00)[1],seq(50,6000,500)),c3=c(coef(a00)[4]/coef(a00)[1],seq(1000,50000,10000))), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0/(1+c1/t+c2*t^(-2) + c3*t^(-3)), data=t2[t2$t>0,], start=coef(a0), control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  #resCoef$sd[i] <- summary(a)$sigma
  tmpHest <- coef(a)[1]/(1 + coef(a)[2]/t2$t + coef(a)[3]/t2$t^2 + coef(a)[4]/t2$t^3)
  tmpHest[!is.finite(tmpHest)] <- 0
  resCoef$sd[i] <- sqrt(sum((t2$h - tmpHest)^2)/(NROW(t2)-4))
   resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2[t2$t>0,], start=coef(a0)[2:4], control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)
      t4 <- t3(ageSteps[j]) / (1/(1+coef(a)[1]/ageSteps[j]+coef(a)[2]*ageSteps[j]^(-2) + coef(a)[3]*ageSteps[j]^(-3)))
      resCoef[i,paste(c("c0","c1","c2","c3"), ageSteps[j], sep="_")] <- c(t4, coef(a))
#      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Hoss1E1"]] <- resCoef
summary(resCoef)


#Yoshida
fn <- function(t, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv / (1/(1+c1/tGiv+c2*tGiv^c3))
  tMul * 1/(1+c1/t+c2*t^c3)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0/(1+c1/t+c2*t^c3), data=t2, start=expand.grid(c0=seq(20,70,10),c1=seq(10,100,10),c2=seq(1000,20000,3000),c3=seq(-3,-1,0.5)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0/(1+c1/t+c2*t^c3), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999, ftol=0.001), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999, ftol=0.001), trace=F)
      t4 <- t3(ageSteps[j]) / (1/(1+coef(a)[1]/ageSteps[j]+coef(a)[2]*ageSteps[j]^coef(a)[3]))
      resCoef[i,paste(c("c0","c1","c2","c3"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Yoshida"]] <- resCoef
summary(resCoef)


#Hossfeld4
#m0 <- nls(H ~ c0/(1+c1*t^c2), start=list(c0=35, c1=8000, c2=-2), data=x)
fn <- function(t, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / (1/(1+c1*tGiv^c2))
  tMul * 1/(1+c1*t^c2)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0/(1+c1*t^c2), data=t2, start=expand.grid(c0=seq(20,80,10),c1=seq(1000,5000,1000),c2=seq(-3,-1,0.25)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0/(1+c1*t^c2), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / (1/(1+coef(a)[1]*ageSteps[j]^coef(a)[2]))
      resCoef[i,paste(c("c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Hossfeld4"]] <- resCoef
summary(resCoef)


#Strand
fn <- function(t, c1, tGiv, hGiv) {
  tMul <- hGiv / (1/(1+c1/tGiv)^3)
  tMul * 1/(1+c1/t)^3
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=6*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1", "hEst", "si"), rep(ageSteps,each=6), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0/(1+c1/t)^3, data=t2, start=expand.grid(c0=seq(20,80,10),c1=seq(1000,5000,1000)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0/(1+c1/t)^3, data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,tGiv,hGiv), data=t2, start=coef(a0)[2], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / (1/(1+coef(a)[1]/ageSteps[j])^3)
      resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Strand"]] <- resCoef
summary(resCoef)



#Levakovic1
#m0 <- nls(H ~ c0/(1+c1*t^-1)^c2, start=list(c0=35, c1=1000, c2=1.5), data=x)
fn <- function(t, c1, c3, tGiv, hGiv) {
  tMul <- hGiv / (1/(1+c1*tGiv^-1)^c3)
  tMul * 1/(1+c1*t^-1)^c3
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c3", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  #a0 <- nls2(h ~ c0/(1+c1*t^-1)^c3, data=t2, start=expand.grid(c0=seq(20,40,5),c1=seq(3000,10000,1000),c3=seq(0.5,1.5,0.2)), trace=F, algorithm = "brute-force")
  a0 <- nls2(h ~ c0/(1+c1*t^-1)^c3, data=t2, start=expand.grid(c0=seq(30,70,10),c1=seq(5,80,10),c3=seq(1,7,1)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0/(1+c1*t^-1)^c3, data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / (1/(1+coef(a)[1]*ageSteps[j]^-1)^coef(a)[2])
      resCoef[i,paste(c("c0","c1","c3"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Levakovic88"]] <- resCoef
summary(resCoef)


#Levakovic2
#m0 <- nls(H ~ c0/(1+c1*t^c3)^c2, start=list(c0=35, c1=1000, c2=1.5,c3=-2), data=x)
fn <- function(t, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv / (1/(1+c1*tGiv^c2)^c3)
  tMul * 1/(1+c1*t^c2)^c3
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0/(1+c1*t^c3)^c2, data=t2, start=expand.grid(c0=seq(30,60,10),c1=c(10,100,300,700,1500,3000,8000,20000),c2=seq(-2.5,-0.5,0.5),c3=seq(0.5,3,0.5)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0/(1+c1*t^c2)^c3, data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / (1/(1+coef(a)[1]*ageSteps[j]^coef(a)[2])^coef(a)[3])
      resCoef[i,paste(c("c0","c1","c2","c3"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Levakovic159"]] <- resCoef
summary(resCoef)


#Logit  -> geht nicht durch 0!!!
#m0 <- nls(H ~ c0/(1+c1*exp(c2*t)), start=list(c0=35, c1=15, c2=-0.05), data=x)
fn <- function(t, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / (1/(1+c1*exp(c2*tGiv)))
  tMul * 1/(1+c1*exp(c2*t))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0/(1+c1*exp(c2*t)), data=t2, start=expand.grid(c0=seq(10,40,10),c1=seq(10,20,2.5),c2=seq(-0.07,-0.03,0.005)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0/(1+c1*exp(c2*t)), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / (1/(1+coef(a)[1]*exp(coef(a)[2]*ageSteps[j])))
      resCoef[i,paste(c("c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Logit"]] <- resCoef
summary(resCoef)


#LogitE1
fn <- function(t, c0, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv - c0/(1+c1*exp(c2*tGiv^c3))
  tMul + c0/(1+c1*exp(c2*t^c3))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric(), c4=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=9*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3","c4", "hEst", "si"), rep(ageSteps,each=9), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nlsLM(h ~ c4 + c0/(1+c1*exp(c2*t)), data=t2, start=list(c0=50, c1=3, c2=-0.03, c4=-10), control=nls.lm.control(maxiter=999), lower=c(-Inf,0,-Inf,-Inf), trace=F)
  a <- nlsLM(h ~ c4 + c0/(1+c1*exp(c2*t^c3)), data=t2, start=c(coef(a0), list(c3=1)), control=nls.lm.control(maxiter=999), lower=c(-Inf,0,-Inf,-Inf,-Inf), trace=F)
  a0 <- a
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c4","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c0,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[c(1:3,5)], control=nls.lm.control(maxiter=999), trace=F)
#      t4 <- t3(ageSteps[j]) / (1/(1+coef(a)[1]*exp(coef(a)[2]*ageSteps[j])))
#      resCoef[i,paste(c("c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["LogitE1"]] <- resCoef
summary(resCoef)


#LogitE2
fn <- function(t, c0, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv - c0/(1+c1*exp(c2*tGiv+c3*tGiv^2))
  tMul + c0/(1+c1*exp(c2*t+c3*t^2))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric(), c4=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=9*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3","c4", "hEst", "si"), rep(ageSteps,each=9), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nlsLM(h ~ c4 + c0/(1+c1*exp(c2*t)), data=t2, start=list(c0=50, c1=3, c2=-0.03, c4=-10), control=nls.lm.control(maxiter=999), lower=c(-Inf,0,-Inf,-Inf), trace=F)
  a <- nlsLM(h ~ c4 + c0/(1+c1*exp(c2*t+c3*t^2)), data=t2, start=c(coef(a0), list(c3=0)), control=nls.lm.control(maxiter=999), lower=c(-Inf,0,-Inf,-Inf,-Inf), trace=F)
  a0 <- a
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c4","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c0,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[c(1:3,5)], control=nls.lm.control(maxiter=999), trace=F)
#      t4 <- t3(ageSteps[j]) / (1/(1+coef(a)[1]*exp(coef(a)[2]*ageSteps[j])))
#      resCoef[i,paste(c("c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["LogitE2"]] <- resCoef
summary(resCoef)



#LogitE3
fn <- function(t, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv / (1/(1+c1*exp(c2*log(tGiv)+c3*log(tGiv)^2)))
  tMul * 1/(1+c1*exp(c2*log(t)+c3*log(t)^2))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0/(1+c1*exp(c2*log(t)+c3*log(t)^2)), data=t2, start=expand.grid(c0=seq(20,70,10),c1=c(200,2000,9000),c2=seq(-4,-1,0.5),c3=0.1), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0/(1+c1*exp(c2*log(t)+c3*log(t)^2)), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999, ftol=0.001), trace=F)
  a0 <- a
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999, ftol=0.001), trace=F)
#      t4 <- t3(ageSteps[j]) / (1/(1+coef(a)[1]*exp(coef(a)[2]*ageSteps[j])))
#      resCoef[i,paste(c("c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["LogitE3"]] <- resCoef
summary(resCoef)


#Nelder
fn <- function(t, c0, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv - c0/(1+c1*exp(c2*tGiv))^c3
  tMul + c0/(1+c1*exp(c2*t))^c3
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric(), c4=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=9*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3","c4", "hEst", "si"), rep(ageSteps,each=9), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nlsLM(h ~ c4 + c0/(1+c1*exp(c2*t)), data=t2, start=list(c0=50, c1=3, c2=-0.03, c4=-10), control=nls.lm.control(maxiter=999), lower=c(-Inf,0,-Inf,-Inf), trace=F)
  a <- tryCatch({nlsLM(h ~ c4 + c0/(1+c1*exp(c2*t))^c3, data=t2, start=c(coef(a0), list(c3=1)), control=nls.lm.control(maxiter=999), lower=c(-Inf,0,-Inf,-Inf,-Inf), trace=F)}, 
                error =function(cond) {nlsLM(h ~ c4 + c0/(1+c1*exp(c2*t))^c3, data=t2, start=c(coef(a0), list(c3=1)), control=nls.lm.control(maxiter=999, ftol=0.0001), lower=c(-Inf,0,-Inf,-Inf,-Inf), trace=F)})
  a0 <- a
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c4","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c0,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[c(1:3,5)], control=nls.lm.control(maxiter=999), trace=F)
#      t4 <- t3(ageSteps[j]) / (1/(1+coef(a)[1]*exp(coef(a)[2]*ageSteps[j])))
#      resCoef[i,paste(c("c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Nelder"]] <- resCoef
summary(resCoef)


#Peschel
fn <- function(t, c1, tGiv, hGiv) {
  tMul <- hGiv / (1 - exp(-2*tGiv/c1) * (1 + 2*tGiv/c1 + 2*(tGiv/c1)^2))
  tMul * (1 - exp(-2*t/c1) * (1 + 2*t/c1 + 2*(t/c1)^2))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=6*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1", "hEst", "si"), rep(ageSteps,each=6), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0 * (1 - exp(-2*t/c1) * (1 + 2*t/c1 + 2*(t/c1)^2)), data=t2, start=expand.grid(c0=seq(10,50,10),c1=seq(20,100,10)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0 * (1 - exp(-2*t/c1) * (1 + 2*t/c1 + 2*(t/c1)^2)), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,tGiv,hGiv), data=t2, start=coef(a0)[2], control=nls.lm.control(maxiter=999), trace=F)
#      t4 <- t3(ageSteps[j]) / (1 - 1/(1+exp(coef(a)[1]*ageSteps[j])))
#      resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Peschel"]] <- resCoef
summary(resCoef)


#Tanh
fn <- function(t, c1,c2,c3, tGiv, hGiv) {
  tMul <- hGiv / tanh(c2*tGiv^c3)^c1
  tMul * tanh(c2*t^c3)^c1
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nlsLM(h ~ c0*tanh(c2*t^c3)^c1, data=t2, start=list(c0=45, c1=10, c2=0.5, c3=0.3), control=nls.lm.control(maxiter=999), trace=F)
  a <- a0
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- tryCatch({nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999), trace=F)},
                    error =function(cond) {a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)})
#      t4 <- t3(ageSteps[j]) / (1 - 1/(1+exp(coef(a)[1]*ageSteps[j])))
#      resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Tanh"]] <- resCoef
summary(resCoef)


#Atan
fn <- function(t, c1,c2,c3, tGiv, hGiv) {
  tMul <- hGiv / atan(c2*tGiv^c3)^c1
  tMul * atan(c2*t^c3)^c1
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nlsLM(h ~ c0*atan(c2*t^c3)^c1, data=t2, start=list(c0=20, c1=2, c2=0.025, c3=1), control=nls.lm.control(maxiter=999), trace=F)
  a <- a0
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- tryCatch({nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999), trace=F)},
                    error =function(cond) {a <- nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999, ftol=0.3), trace=F)})
#      t4 <- t3(ageSteps[j]) / (1 - 1/(1+exp(coef(a)[1]*ageSteps[j])))
#      resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Atan"]] <- resCoef
summary(resCoef)


#Log
fn <- function(t, c1,c2,c3, tGiv, hGiv) {
  tMul <- hGiv / log(1 + exp(c2)*tGiv^c3)^c1
  tMul * log(1 + exp(c2)*t^c3)^c1
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^c3)^c1, data=t2, start=list(c0=15, c1=0.6, c2=-10, c3=3), control=nls.lm.control(maxiter=999), trace=F)
  a <- a0
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- tryCatch({nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999), lower=c(0,0,0), trace=F)},
                    error =function(cond) {nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(a0)[2:4], control=nls.lm.control(maxiter=999, ftol=0.001), lower=c(0,0,0), trace=F)})
#      t4 <- t3(ageSteps[j]) / (1 - 1/(1+exp(coef(a)[1]*ageSteps[j])))
#      resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Ln"]] <- resCoef
summary(resCoef)


#Koevessi1
fn <- function(t, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv / ((1 - exp(c1*tGiv)) + c2*(1 - exp(c3*tGiv)))
  tMul * ((1 - exp(c1*t)) + c2*(1 - exp(c3*t)))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0*(1 - exp(c1*t)) + c2*(1 - exp(c3*t)), data=t2, start=expand.grid(c0=seq(30,200,20),c1=seq(-0.04,-0.01,0.01),c2=seq(-100,-10,20),c3=seq(-0.2,-0.02,0.04)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0*(1 - exp(c1*t)) + c2*(1 - exp(c3*t)), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  a00 <- a
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      tCoef <- coef(a00)[2:4]; tCoef[2] <- tCoef[2] / coef(a00)[1]
      tCoef0 <- coef(a0)[2:4]; tCoef0[2] <- tCoef0[2] / coef(a0)[1]
      a <- tryCatch({nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=tCoef, control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)},
                    error =function(cond) {nlsLM(h ~ fn(t,c1,c2,c3,tGiv,hGiv), data=t2, start=tCoef0, control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)})
      t4 <- t3(ageSteps[j]) / ((1 - exp(coef(a)[1]*ageSteps[j])) + coef(a)[2]*(1 - exp(coef(a)[3]*ageSteps[j])))
      resCoef[i,paste(c("c0","c1","c2","c3"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Koevessi1"]] <- resCoef
summary(resCoef)

#Koevessi2
fn <- function(t, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / ((1-exp(-c1*tGiv))/c1 - (1-exp(-c2*tGiv))/c2)
  tMul * ((1-exp(-c1*t))/c1 - (1-exp(-c2*t))/c2)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a0 <- nls2(h ~ c0 * ((1-exp(-c1*t))/c1 - (1-exp(-c2*t))/c2), data=t2, start=expand.grid(c0=seq(1,7,1),c1=seq(0.01,0.04,0.01),c2=seq(0.015,0.2,0.01)), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ c0 * ((1-exp(-c1*t))/c1 - (1-exp(-c2*t))/c2), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), lower=c(0,0,0), upper=c(100,1,1), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.lm.control(maxiter=999, ftol=0.001), trace=F)
      t4 <- t3(ageSteps[j]) / ((1-exp(-coef(a)[1]*ageSteps[j]))/coef(a)[1] - (1-exp(-coef(a)[2]*ageSteps[j]))/coef(a)[2])
      resCoef[i,paste(c("c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Koevessi2"]] <- resCoef
summary(resCoef)


##Grosenbaugh
fn <- function(t, c1, c2, c3, c4, tGiv, hGiv) {
  tMul <- hGiv / ((exp((c1^2-1)*exp(c2*tGiv^c3)) - c1 * exp(c2*tGiv^c3))^(c1*c4+1))
  tMul * ((exp((c1^2-1)*exp(c2*t^c3)) - c1 * exp(c2*t^c3))^(c1*c4+1))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric(), c4=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=9*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3","c4", "hEst", "si"), rep(ageSteps,each=9), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  a <- tryCatch({nlsLM(h ~ c0*(exp((c1^2-1)*exp(c2*t^c3)) - c1 * exp(c2*t^c3))^(c1*c4+1), data=t2, start=list(c0=45, c1=-9, c2=-1, c3=0.3, c4=0.2), control=nls.lm.control(maxiter=999), trace=F)},
                error =function(cond) {nlsLM(h ~ c0*(exp((c1^2-1)*exp(c2*t^c3)) - c1 * exp(c2*t^c3))^(c1*c4+1), data=t2, start=list(c0=45, c1=-9, c2=-1, c3=0.3, c4=0.2), control=nls.lm.control(maxiter=999, ftol=0.1), trace=F)})
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3","c4")] <- coef(a)
  a00 <- a
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- a <- tryCatch({nlsLM(h ~ fn(t,c1,c2,c3,c4,tGiv,hGiv), data=t2, start=coef(a00)[2:5], control=nls.lm.control(maxiter=999), trace=F)},
                         error =function(cond) {nls(h ~ fn(t,c1,c2,c3,c4,tGiv,hGiv), data=t2, start=coef(a00)[2:5], control=nls.control(maxiter=150, warnOnly=T), algo="port", trace=F)})
#      t4 <- t3(ageSteps[j]) / (1 - 1/(1+exp(coef(a)[1]*ageSteps[j])))
#      resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
    }
  }
}
res[["Grosenbaugh"]] <- resCoef
summary(resCoef)



##Backman
f1 <- function(h, dt, c0, c1, c2) {
  h <- rep(h,length(dt))
  for(i in 1:length(h)) {
    if(dt[i] > 0) {
      for(j in 1:dt[i]) {
        h[i] <- h[i] + c0*exp(c1*log(j) + c2*log(j)^2)
      }
    }
  }
  h
}
fn <- function(h, t, c0, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / f1(h, tGiv, c0, c1, c2)
  tMul * f1(h, t, c0, c1, c2)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c00","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
i <- 100
for(i in 1:length(t1)) {
  print(i)
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  t2 <- t2[order(t2$t),]
  t3 <- data.frame(dt = diff(t2$t), ih = diff(t2$h), t=tail(t2$t, -1))
  #Find good starting parameters 
  a00 <- lm(log(ih/dt) ~ log(t) + I(log(t)^2), data=t3)
  a0 <- nls2(ih/dt ~ c0*exp(c1*log(t) + c2*log(t)^2), data=t3, start=expand.grid(c0=exp(coef(a00)[1]),c1=coef(a00)[2],c2=coef(a00)[3]), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ f1(0, t, c0, c1, c2), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  t5 <- a
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(0, t,c0,c1,c2,tGiv,hGiv), data=t2, start=coef(t5)[1:3], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / f1(0, ageSteps[j], coef(a)[1], coef(a)[2], coef(a)[3])
      resCoef[i,paste(c("c00", "c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100, tGiv=ageSteps[j], hGiv=t3(ageSteps[j])))
    }
  }
}
res[["Backman"]] <- resCoef
summary(resCoef)


##BackmanE1
f1 <- function(h, dt, c0, c1, c2, c3) {
  h <- rep(h,length(dt))
  for(i in 1:length(h)) {
    if(dt[i] > 0) {
      for(j in 1:dt[i]) {
        h[i] <- h[i] + c0*exp(c1*log(j) + c2*log(j)^2 + c3*log(j)^4)
      }
    }
  }
  h
}
fn <- function(h, t, c0, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv / f1(h, tGiv, c0, c1, c2, c3)
  tMul * f1(h, t, c0, c1, c2, c3)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=9*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c00","c0","c1","c2","c3","hEst", "si"), rep(ageSteps,each=9), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
i <- 63
for(i in 1:length(t1)) {
  print(i)
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  t2 <- t2[order(t2$t),]
  t3 <- data.frame(dt = diff(t2$t), ih = diff(t2$h), t=tail(t2$t, -1))
  #Find good starting parameters 
  a00 <- lm(log(ih/dt) ~ log(t) + I(log(t)^2) + I(log(t)^4), data=t3)
  #a0 <- nls2(ih/dt ~ c0*exp(c1*log(t) + c2*log(t)^2 + c3*log(t)^4), data=t3, start=expand.grid(c0=c(10, exp(coef(a00)[1])),c1=c(-10, coef(a00)[2]),c2=c(0.5, coef(a00)[3]),c3=c(-0.05, coef(a00)[4])), trace=F, algorithm = "brute-force")
  a0 <- nls2(h ~ f1(0, t, c0, c1, c2, c3), data=t2, start=expand.grid(c0=c(0.00002,10, exp(coef(a00)[1])),c1=c(7,-10, coef(a00)[2]),c2=c(-1,0.5, coef(a00)[3]),c3=c(0.005,-0.05, coef(a00)[4])), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ f1(0, t, c0, c1, c2, c3), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999, ftol=0.001), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  t5 <- a
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(0, t,c0,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(t5)[1:4], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / f1(0, ageSteps[j], coef(a)[1], coef(a)[2], coef(a)[3], coef(a)[4])
      resCoef[i,paste(c("c00", "c0","c1","c2","c3"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100, tGiv=ageSteps[j], hGiv=t3(ageSteps[j])))
    }
  }
}
res[["BackmanE1"]] <- resCoef
summary(resCoef)


##Koller
f1 <- function(h, dt, c0, c1, c2) {
  h <- rep(h,length(dt))
  for(i in 1:length(h)) {
    if(dt[i] > 0) {
      for(j in 1:dt[i]) {
        h[i] <- h[i] + c0*j^c1*c2^(-j)
      }
    }
  }
  h
}
fn <- function(h, t, c0, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / f1(h, tGiv, c0, c1, c2)
  tMul * f1(h, t, c0, c1, c2)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c00","c0","c1","c2", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
i <- 100
for(i in 1:length(t1)) {
  print(i)
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  t2 <- t2[order(t2$t),]
  t3 <- data.frame(dt = diff(t2$t), ih = diff(t2$h), t=tail(t2$t, -1))
  #Find good starting parameters 
  a00 <- lm(log(ih/dt) ~ log(t) + I(-t), data=t3)
  a0 <- nls2(ih/dt ~ c0*t^c1*c2^(-t), data=t3, start=expand.grid(c0=exp(coef(a00)[1]),c1=coef(a00)[2],c2=exp(coef(a00)[3])), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ f1(0, t, c0, c1, c2), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  t5 <- a
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(0, t,c0,c1,c2,tGiv,hGiv), data=t2, start=coef(t5)[1:3], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / f1(0, ageSteps[j], coef(a)[1], coef(a)[2], coef(a)[3])
      resCoef[i,paste(c("c00", "c0","c1","c2"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100, tGiv=ageSteps[j], hGiv=t3(ageSteps[j])))
    }
  }
}
res[["Koller"]] <- resCoef
summary(resCoef)


##Levakovic147
f1 <- function(h, dt, c0, c1, c2, c3) {
  h <- rep(h,length(dt))
  for(i in 1:length(h)) {
    if(dt[i] > 0) {
      for(j in 1:dt[i]) {
        h[i] <- h[i] + c0*j^c1*c2^(-j^c3)
      }
    }
  }
  h
}
fn <- function(h, t, c0, c1, c2, c3, tGiv, hGiv) {
  tMul <- hGiv / f1(h, tGiv, c0, c1, c2, c3)
  tMul * f1(h, t, c0, c1, c2, c3)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=9*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c00","c0","c1","c2","c3","hEst", "si"), rep(ageSteps,each=9), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
i <- 100
for(i in 1:length(t1)) {
  print(i)
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  t2 <- t2[order(t2$t),]
  t3 <- data.frame(dt = diff(t2$t), ih = diff(t2$h), t=tail(t2$t, -1))
  #Find good starting parameters 
  a00 <- lm(log(ih/dt) ~ log(t) + I(-t), data=t3)
  a0 <- nls2(ih/dt ~ c0*t^c1*c2^(-t), data=t3, start=expand.grid(c0=exp(coef(a00)[1]),c1=coef(a00)[2],c2=exp(coef(a00)[3]), c3=1), trace=F, algorithm = "brute-force")
  a <- nlsLM(h ~ f1(0, t, c0, c1, c2,c3), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  t5 <- a
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- nlsLM(h ~ fn(0, t,c0,c1,c2,c3,tGiv,hGiv), data=t2, start=coef(t5)[1:4], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / f1(0, ageSteps[j], coef(a)[1], coef(a)[2], coef(a)[3], coef(a)[4])
      resCoef[i,paste(c("c00", "c0","c1","c2","c3"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100, tGiv=ageSteps[j], hGiv=t3(ageSteps[j])))
    }
  }
}
res[["Levakovic147"]] <- resCoef
summary(resCoef)


##Hyperlogistic
f1 <- function(h, dt, c0, c1, c2, c3, c4) {
  h <- rep(h,length(dt))
  for(i in 1:length(h)) {
    if(c2 > (c4+h[i])) {
      if(dt[i] > 0) {
        for(j in 1:dt[i]) {
          h[i] <- min(h[i] + (c0*(c4+h[i])^c1)*(1-(c4+h[i])/c2)^c3, c2)
        }
      }
    }
  }
  h
}
fn <- function(h, t, c0, c1, c2, c3, c4, tGiv, hGiv) {
  tMul <- hGiv / f1(h, tGiv, c0, c1, c2, c3, c4)
  tMul * f1(h, t, c0, c1, c2, c3, c4)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric(), c4=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=10*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c00","c0","c1","c2","c3","c4", "hEst", "si"), rep(ageSteps,each=10), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
i <- 4
for(i in 1:length(t1)) {
  print(i)
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  t2 <- t2[order(t2$t),]
  t3 <- data.frame(dt = diff(t2$t), ih = diff(t2$h), h=head(t2$h, -1))
  #Find good starting parameters 
  a0 <- nlsLM(ih/dt ~ (c0*(c4+h)^c1)*(1-(c4+h)/c2)^c3, data=t3, start=list(c0=0.2, c1=0.4, c2=500, c3=35, c4=1.5), control=nls.lm.control(maxiter=999), lower=c(0.05,0.1,max(t2$h)*1.3,0.01,0), upper=c(1,5,Inf,Inf,5), trace=F)
  a <- tryCatch({nlsLM(h ~ f1(0, t, c0, c1, c2, c3, c4), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), lower=c(0,0,0,0,0), trace=F)},
                error =function(cond) {nlsLM(h ~ f1(0, t, c0, c1, c2, c3, c4), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999, ftol=0.001), lower=c(0,0.01,max(t2$h)+2,0.01,0), trace=F)})
  a1 <- tryCatch({nlsLM(h ~ f1(0, t, c0, c1, c2, c3, c4), data=t2, start=list(c0=0.2, c1=0.4, c2=500, c3=35, c4=1.5), control=nls.lm.control(maxiter=999, ftol=0.001), lower=c(0,0,0,0,0), trace=F)},
                 error =function(cond) {nlsLM(h ~ f1(0, t, c0, c1, c2, c3, c4), data=t2, start=list(c0=0.3, c1=0.3, c2=100, c3=30, c4=0.001), control=nls.lm.control(maxiter=999, ftol=0.001), lower=c(0.05,0.1,max(t2$h)*1.3,0.01,0), trace=F)})
  a2 <- nlsLM(h ~ f1(0, t, c0, c1, c2, c3, c4), data=t2, start=list(c0=0.3, c1=0.3, c2=100, c3=30, c4=0.001), control=nls.lm.control(maxiter=999, ftol=0.001), lower=c(0.05,0.1,max(t2$h)*1.3,0.01,0), trace=F)
  if(summary(a)$sigma > summary(a1)$sigma) {a <- a1}
  if(summary(a)$sigma > summary(a2)$sigma) {a <- a2}
  a0 <- a
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3","c4")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  t5 <- a
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- tryCatch({nlsLM(h ~ fn(0, t,coef(t5)[1],c1,c2,c3,c4,tGiv,hGiv), data=t2, start=coef(t5)[2:5], control=nls.lm.control(maxiter=999), trace=F)},
                    error =function(cond) {nlsLM(h ~ fn(0, t,coef(t5)[1],c1,c2,c3,c4,tGiv,hGiv), data=t2, start=coef(t5)[2:5], control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)})
#      t4 <- t3(ageSteps[j]) / f1(0, ageSteps[j], coef(t5)[1], coef(a)[2], coef(a)[3], coef(a)[4], coef(a)[5])
#      resCoef[i,paste(c("c00", "c0","c1","c2","c3","c4"), ageSteps[j], sep="_")] <- c(t4, coef(t5)[1], coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100, tGiv=ageSteps[j], hGiv=t3(ageSteps[j])))
    }
  }
}
res[["HyperLog"]] <- resCoef
summary(resCoef)


##EvolonE1
f1 <- function(h, dt, c0, c1, c2, c3, c4) {
  h <- rep(h,length(dt))
  for(i in 1:length(h)) {
    if(c2 > h[i]) {
      if(dt[i] > 0) {
        for(j in 1:dt[i]) {
          h[i] <- min(h[i] + (c4 + c0*h[i]^c1)*(1-h[i]/c2)^c3, c2)
        }
      }
    }
  }
  h
}
fn <- function(h, t, c0, c1, c2, c3, c4, tGiv, hGiv) {
  tMul <- hGiv / f1(h, tGiv, c0, c1, c2, c3, c4)
  tMul * f1(h, t, c0, c1, c2, c3, c4)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric(), c4=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=10*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c00","c0","c1","c2","c3","c4", "hEst", "si"), rep(ageSteps,each=10), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
i <- 147
i <- 56
for(i in 1:length(t1)) {
  print(i)
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  t2 <- t2[order(t2$t),]
  t3 <- data.frame(dt = diff(t2$t), ih = diff(t2$h), h=head(t2$h, -1))
  #Find good starting parameters 
  a0 <- nlsLM(ih/dt ~ (c4+c0*h^c1)*(1-h/c2)^c3, data=t3, start=list(c0=0.15, c1=0.65, c2=200, c3=1, c4=0.10), control=nls.lm.control(maxiter=999), lower=c(0.05,0.1,10,0.01,0), upper=c(1,5,Inf,Inf,1), trace=F)
  if(i==147) {
    a <- nlsLM(h ~ f1(0, t, c0, c1, c2, c3, c4), data=t2, start=list(c0=0.25, c1=0.37, c2=77, c3=3.5, c4=8e-312), control=nls.lm.control(maxiter=999), trace=F)
  } else {
    a <- nlsLM(h ~ f1(0, t, c0, c1, c2, c3, c4), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  }
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3","c4")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  t5 <- a
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      if(i %in% c(56,125,147)) {
        a <- nlsLM(h ~ fn(0, t,coef(t5)[1],c1,c2,c3,c4,tGiv,hGiv), data=t2, start=list(c1=.4,c2=80,c3=3.5,c4=0.01), control=nls.lm.control(maxiter=999), trace=F)
      } else {
        a <- nlsLM(h ~ fn(0, t,coef(t5)[1],c1,c2,c3,c4,tGiv,hGiv), data=t2, start=coef(t5)[2:5], control=nls.lm.control(maxiter=999), trace=F)
      }
      t4 <- t3(ageSteps[j]) / f1(0, ageSteps[j], coef(t5)[1], coef(a)[2], coef(a)[3], coef(a)[4], coef(a)[5])
      resCoef[i,paste(c("c00", "c0","c1","c2","c3","c4"), ageSteps[j], sep="_")] <- c(t4, coef(t5)[1], coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100, tGiv=ageSteps[j], hGiv=t3(ageSteps[j])))
    }
  }
}
res[["HyperLogE1"]] <- resCoef
summary(resCoef)


##Kinder
#mF <- nls(ih/U ~ c0*(c1+H)^c2*exp(c3*H^c4), data=baum, start=list(c0=3, c1=0.01, c2=0.3, c3=-0.01, c4=1), algo="port", trace=T, control=nls.control(maxit=999, warnOnly = T), weight=baum$U)
f1 <- function(h, dt, c0, c1, c2, c3, c4) {
  h <- rep(h,length(dt))
  for(i in 1:length(h)) {
    if(dt[i] > 0) {
      for(j in 1:dt[i]) {
        h[i] <- h[i] + (c0 + c1*h[i]^c2)*exp(c3*h[i]^c4)
      }
    }
  }
  h
}
fn <- function(h, t, c0, c1, c2, c3, c4, tGiv, hGiv) {
  tMul <- hGiv / f1(h, tGiv, c0, c1, c2, c3, c4)
  tMul * f1(h, t, c0, c1, c2, c3, c4)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric(), c4=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=10*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c00","c0","c1","c2","c3","c4", "hEst", "si"), rep(ageSteps,each=10), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
i <- 82
for(i in 1:length(t1)) {
  print(i)
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  t2 <- t2[order(t2$t),]
  t3 <- data.frame(dt = diff(t2$t), ih = diff(t2$h), h=head(t2$h, -1))
  #Find good starting parameters 
  a0 <- nlsLM(ih/dt ~ (c0+c1*h^c2)*exp(c3*h^c4), data=t3, start=expand.grid(c0=0.15,c1=0.25,c2=1, c3=-0.3, c4=0.7), control=nls.lm.control(maxiter=999, ftol=0.01), trace=F)
  a <- nlsLM(h ~ f1(0, t, c0, c1, c2, c3, c4), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  #plot(t2$t, t2$h, type="b")
  #lines(t2$t, predict(a), col=2)
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3","c4")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  t5 <- a
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      #a <- nlsLM(h ~ fn(0, t,c0,c1,c2,c3,c4,tGiv,hGiv), data=t2, start=coef(t5)[1:5], control=nls.lm.control(maxiter=999), trace=F)
      #t4 <- t3(ageSteps[j]) / f1(0, ageSteps[j], coef(a)[1], coef(a)[2], coef(a)[3], coef(a)[4], coef(a)[5])
      #resCoef[i,paste(c("c00", "c0","c1","c2","c3","c4"), ageSteps[j], sep="_")] <- c(t4, coef(a))
      a <- nlsLM(h ~ fn(0, t,coef(t5)[1],c1,c2,c3,c4,tGiv,hGiv), data=t2, start=coef(t5)[2:5], control=nls.lm.control(maxiter=999), trace=F)
      t4 <- t3(ageSteps[j]) / f1(0, ageSteps[j], coef(t5)[1], coef(a)[2], coef(a)[3], coef(a)[4], coef(a)[5])
      resCoef[i,paste(c("c00", "c0","c1","c2","c3","c4"), ageSteps[j], sep="_")] <- c(t4, coef(t5)[1], coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100, tGiv=ageSteps[j], hGiv=t3(ageSteps[j])))
    }
  }
}
res[["HyperLogE2"]] <- resCoef
summary(resCoef)


##HyperlogE3
f1 <- function(h, dt, c0, c1, c2, c3, c4) {
  h <- rep(h,length(dt))
  for(i in 1:length(h)) {
    if(dt[i] > 0) {
      for(j in 1:dt[i]) {
        h[i] <- h[i] + (c0 + c1*h[i]^c2)*(1/(1 + c3*h[i]^c4))
      }
    }
  }
  h
}
fn <- function(h, t, c0, c1, c2, c3, c4, tGiv, hGiv) {
  tMul <- hGiv / f1(h, tGiv, c0, c1, c2, c3, c4)
  tMul * f1(h, t, c0, c1, c2, c3, c4)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric(), c4=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=10*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c00","c0","c1","c2","c3","c4", "hEst", "si"), rep(ageSteps,each=10), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
i <- 147
for(i in 1:length(t1)) {
  print(i)
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  t2 <- t2[order(t2$t),]
  t3 <- data.frame(dt = diff(t2$t), ih = diff(t2$h), h=head(t2$h, -1))
  #Find good starting parameters 
  a0 <- tryCatch({nlsLM(ih/dt ~ (c0+c1*h^c2)*(1/(1 + c3*h^c4)), data=t3, start=list(c0=0.02, c1=0.15, c2=0.5, c3=0.0061, c4=3), control=nls.lm.control(maxiter=999, ftol=0.1), lower=c(1e-5,1e-5,1e-5,1e-5,1e-5), trace=F)},
                 error =function(cond) {nlsLM(ih/dt ~ (c0+c1*h^c2)*(1/(1 + c3*h^c4)), data=t3, start=list(c0=0.2, c1=0.15, c2=0.5, c3=0.0061, c4=3), control=nls.lm.control(maxiter=999, ftol=0.1), lower=c(1e-5,1e-5,1e-5,1e-5,1e-5), trace=F)})
  a <- tryCatch({nlsLM(h ~ f1(0, t, c0, c1, c2, c3, c4), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), lower=c(1e-5,1e-5,1e-5,1e-5,1e-5), trace=F)},
                error =function(cond) {nlsLM(h ~ f1(0, t, c0, c1, c2, c3, c4), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999, ftol=0.01), lower=c(1e-5,1e-5,1e-5,1e-5,1e-5), trace=F)})
  a0 <- a
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3","c4")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  t5 <- a
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- tryCatch({nlsLM(h ~ fn(0, t,coef(t5)[1],c1,c2,c3,c4,tGiv,hGiv), data=t2, start=coef(t5)[2:5], control=nls.lm.control(maxiter=999, ftol=0.1), trace=F)},
                    error=function(cond) {nlsLM(h ~ fn(0, t,c0,c1,c2,c3,c4,tGiv,hGiv), data=t2, start=coef(t5), control=nls.lm.control(maxiter=999, ftol=0.1), trace=F)})
#      t4 <- t3(ageSteps[j]) / f1(0, ageSteps[j], coef(t5)[1], coef(a)[2], coef(a)[3], coef(a)[4], coef(a)[5])
#      resCoef[i,paste(c("c00", "c0","c1","c2","c3","c4"), ageSteps[j], sep="_")] <- c(t4, coef(t5)[1], coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100, tGiv=ageSteps[j], hGiv=t3(ageSteps[j])))
    }
  }
}
res[["HyperLogE3"]] <- resCoef
summary(resCoef)


##HyperlogE4
f1 <- function(h, dt, c0, c1, c2, c3, c4) {
  h <- rep(h,length(dt))
  for(i in 1:length(h)) {
    if(dt[i] > 0) {
      for(j in 1:dt[i]) {
        h[i] <- h[i] + (c0 + c1*h[i]^c2)*(1 - tanh(c3*h[i]^c4))
      }
    }
  }
  h
}
fn <- function(h, t, c0, c1, c2, c3, c4, tGiv, hGiv) {
  tMul <- hGiv / f1(h, tGiv, c0, c1, c2, c3, c4)
  tMul * f1(h, t, c0, c1, c2, c3, c4)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric(), c4=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=10*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c00","c0","c1","c2","c3","c4", "hEst", "si"), rep(ageSteps,each=10), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
i <- 147
for(i in 1:length(t1)) {
  print(i)
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  t2 <- t2[order(t2$t),]
  t3 <- data.frame(dt = diff(t2$t), ih = diff(t2$h), h=head(t2$h, -1))
  #Find good starting parameters 
  a0 <- nlsLM(ih/dt ~ (c0+c1*h^c2)*(1 - tanh(c3*h^c4)), data=t3, start=list(c0=0.1, c1=0.2, c2=0.8, c3=0.06, c4=1), control=nls.lm.control(maxiter=999, ftol=0.1), lower=c(1e-5,1e-5,1e-5,1e-5,1e-5), trace=F)
  a <- tryCatch({nlsLM(h ~ f1(0, t, c0, c1, c2, c3, c4), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), lower=c(1e-5,1e-5,1e-5,1e-5,1e-5), trace=F)},
                error =function(cond) {nlsLM(h ~ f1(0, t, c0, c1, c2, c3, c4), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999, ftol=0.01), lower=c(1e-5,1e-5,1e-5,1e-5,1e-5), trace=F)})
  a0 <- a
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3","c4")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  t5 <- a
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- tryCatch({nlsLM(h ~ fn(0, t,coef(t5)[1],c1,c2,c3,c4,tGiv,hGiv), data=t2, start=coef(t5)[2:5], control=nls.lm.control(maxiter=999, ftol=0.1), trace=F)},
                    error=function(cond) {nlsLM(h ~ fn(0, t,c0,c1,c2,c3,c4,tGiv,hGiv), data=t2, start=coef(t5), control=nls.lm.control(maxiter=999, ftol=0.1), trace=F)})
#      t4 <- t3(ageSteps[j]) / f1(0, ageSteps[j], coef(t5)[1], coef(a)[2], coef(a)[3], coef(a)[4], coef(a)[5])
#      resCoef[i,paste(c("c00", "c0","c1","c2","c3","c4"), ageSteps[j], sep="_")] <- c(t4, coef(t5)[1], coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100, tGiv=ageSteps[j], hGiv=t3(ageSteps[j])))
    }
  }
}
res[["HyperLogE4"]] <- resCoef
summary(resCoef)


##HyperlogE5
f1 <- function(h, dt, c0, c1, c2, c3, c4) {
  h <- rep(h,length(dt))
  for(i in 1:length(h)) {
    if(dt[i] > 0) {
      for(j in 1:dt[i]) {
        h[i] <- h[i] + (c0 + c1*h[i]^c2)*(pi/2 - atan(c3*h[i]^c4))
      }
    }
  }
  h
}
fn <- function(h, t, c0, c1, c2, c3, c4, tGiv, hGiv) {
  tMul <- hGiv / f1(h, tGiv, c0, c1, c2, c3, c4)
  tMul * f1(h, t, c0, c1, c2, c3, c4)
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric(), c4=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=10*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c00","c0","c1","c2","c3","c4", "hEst", "si"), rep(ageSteps,each=10), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
i <- 147
for(i in 1:length(t1)) {
  print(i)
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  t2 <- t2[order(t2$t),]
  t3 <- data.frame(dt = diff(t2$t), ih = diff(t2$h), h=head(t2$h, -1))
  #Find good starting parameters 
  a0 <- nlsLM(ih/dt ~ (c0+c1*h^c2)*(pi/2 - atan(c3*h^c4)), data=t3, start=list(c0=0.01, c1=0.1, c2=0.5, c3=0.001, c4=2.5), control=nls.lm.control(maxiter=999, ftol=0.1), lower=c(1e-5,1e-5,1e-5,1e-5,1e-5), trace=F)
  a <- tryCatch({nlsLM(h ~ f1(0, t, c0, c1, c2, c3, c4), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), lower=c(1e-5,1e-5,1e-5,1e-5,1e-5), trace=F)},
                error =function(cond) {nlsLM(h ~ f1(0, t, c0, c1, c2, c3, c4), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999, ftol=0.01), lower=c(1e-5,1e-5,1e-5,1e-5,1e-5), trace=F)})
  a0 <- a
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2","c3","c4")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
  t5 <- a
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 10) { #nicht extrapolieren
      t2$tGiv <- ageSteps[j]
      t2$hGiv <- t3(ageSteps[j])
      a <- tryCatch({nlsLM(h ~ fn(0, t,coef(t5)[1],c1,c2,c3,c4,tGiv,hGiv), data=t2, start=coef(t5)[2:5], control=nls.lm.control(maxiter=999, ftol=0.1), trace=F)},
                    error=function(cond) {nlsLM(h ~ fn(0, t,c0,c1,c2,c3,c4,tGiv,hGiv), data=t2, start=coef(t5), control=nls.lm.control(maxiter=999, ftol=0.1), trace=F)},
                    error=function(cond) {nlsLM(h ~ fn(0, t,0.1,c1,c2,c3,c4,tGiv,hGiv), data=t2, start=coef(t5)[2:5], control=nls.lm.control(maxiter=999, ftol=0.1), trace=F)})
#      t4 <- t3(ageSteps[j]) / f1(0, ageSteps[j], coef(t5)[1], coef(a)[2], coef(a)[3], coef(a)[4], coef(a)[5])
#      resCoef[i,paste(c("c00", "c0","c1","c2","c3","c4"), ageSteps[j], sep="_")] <- c(t4, coef(t5)[1], coef(a))
      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100, tGiv=ageSteps[j], hGiv=t3(ageSteps[j])))
    }
  }
}
res[["HyperLogE5"]] <- resCoef
summary(resCoef)


#Thomasius
fn <- function(t, c1, c2, tGiv, hGiv) {
  tMul <- hGiv / ((1 - exp(c1*tGiv*(1-exp(c2*tGiv)))))
  tMul * (1 - exp(c1*t*(1-exp(c2*t))))
}
t1 <- unique(x$id)  #Baumidentifiers
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=7*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","hEst", "si"), rep(ageSteps,each=7), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
for(i in 1:length(t1)) {
  t2 <- x[x$id == t1[i],]  #Data of single tree
  #t2 <- t2[t2$t>0,] #data h=0, t=0 does not help for this equation
  #Find good starting parameters 
  #a0 <- nls2(h ~ c0*(1 - exp(c1*t*(1-exp(c2*t)))), data=t2, start=expand.grid(c0=seq(10,80,10),c1=seq(-0.035,-0.005,0.005), c2=-1), trace=F, algorithm = "brute-force")
  #a <- nlsLM(h ~ c0*(1 - exp(c1*t*(1-exp(c2*t)))), data=t2, start=coef(a0), control=nls.lm.control(maxiter=999), trace=F)
  a <- nlsLM(h ~ c0*(1 - exp(c1*t*(1-exp(c2*t)))), data=t2, start=c(c0=40, c1=-0.02, c2=-0.07), control=nls.lm.control(maxiter=999), trace=F, upper=c(120,0.0001,0.0001), lower=c(0,-0.2,-1))
  a0 <- 1
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
  resCoef[i,c("c0","c1","c2")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
#  for(j in 1:length(ageSteps)) {
#    if(max(t2$t) >= ageSteps[j] && ageSteps[j] > 0) { #nicht extrapolieren
#      t2$tGiv <- ageSteps[j]
#      t2$hGiv <- t3(ageSteps[j])
#      a <- nlsLM(h ~ fn(t,c1,c2,tGiv,hGiv), data=t2, start=coef(a0)[2:3], control=nls.lm.control(maxiter=999), trace=F)
#      t4 <- t3(ageSteps[j]) / ((1 - exp(coef(a)[1]*ageSteps[j])))
#      resCoef[i,paste(c("c0","c1"), ageSteps[j], sep="_")] <- c(t4, coef(a))
#      resCoef[i,paste("sd", ageSteps[j], sep="_")] <- summary(a)$sigma
#      resCoef[i,paste("si", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=100))[1]
#    }
#  }
}
res[["Thomasius"]] <- resCoef
summary(resCoef)



#order(resCoef$sd)


save(res, file="results.RData")

