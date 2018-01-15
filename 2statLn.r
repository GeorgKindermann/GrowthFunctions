#R <2statLn.r --no-save

library(minpack.lm)

datStam <- read.csv("cleanData.csv")
ageSteps <- c(0, 0.0001, 10, 30, 70, 100, 150, 300, 800)
#res <- list()
load("results.RData")


t1 <- unique(datStam$id)
resCoef <- data.frame(src=character(), id=character(), si=numeric(), sd=numeric(), c0=numeric(), c1=numeric(), c2=numeric(), c3=numeric())
t2 <- as.data.frame(matrix(numeric(), ncol=8*length(ageSteps), nrow=0))
names(t2) <- paste(c("h", "sd","c0","c1","c2","c3", "hEst", "si"), rep(ageSteps,each=8), sep="_")
resCoef <- cbind(resCoef, t2)
resCoef[1:length(t1),] <- NA
resCoef$id <- t1
resCoef$src <- as.character(resCoef$src)
resCoefTempl <- resCoef


#h = c0*log(1 + exp(c2)*t^c3)^c1fix
fun <- function(par, data) {
  c1 <- par
  t1 <- unique(data$id)
  resSd <- rep(NA, length(t1))
  for(i in 1:length(t1)) {
    t2 <- data[data$id == t1[i],]  #Data of single tree
    a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^c3)^c1, data=t2, start=list(c0=15, c2=-10, c3=3), control=nls.lm.control(maxiter=999), trace=F)
    resSd[i] <- summary(a0)$sigma
  }
  median(resSd)
}
cfix <- c(optimize(fun, interval=c(0.4, 0.65), data=datStam[datStam$src=="gut",])$minimum, optimize(fun, interval=c(0.4, 0.65), data=datStam[datStam$src=="jrm",])$minimum)
names(cfix) <- c("gut", "jrm")

resCoef <- resCoefTempl
t1 <- unique(datStam$id)
for(i in 1:length(t1)) {
  t2 <- datStam[datStam$id == t1[i],]  #Data of single tree
  c1 <- cfix[t2$src[1]]
  a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^c3)^c1, data=t2, start=list(c0=15, c2=-10, c3=3), control=nls.lm.control(maxiter=999), trace=F)
  a <- a0
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
#  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
}
res[["LnFixC1"]] <- resCoef
summary(resCoef)

#h = c0*log(1 + exp(c2)*t^c3fix)^c1
fun <- function(par, data) {
  c3 <- par
  t1 <- unique(data$id)
  resSd <- rep(NA, length(t1))
  for(i in 1:length(t1)) {
    t2 <- data[data$id == t1[i],]  #Data of single tree
    a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^c3)^c1, data=t2, start=list(c0=15, c1=0.6, c2=-10), control=nls.lm.control(maxiter=999), trace=F)
    resSd[i] <- summary(a0)$sigma
  }
  median(resSd)
}
cfix <- c(optimize(fun, interval=c(2, 4), data=datStam[datStam$src=="gut",])$minimum, optimize(fun, interval=c(2, 4), data=datStam[datStam$src=="jrm",])$minimum)
names(cfix) <- c("gut", "jrm")
resCoef <- resCoefTempl
t1 <- unique(datStam$id)
for(i in 1:length(t1)) {
  t2 <- datStam[datStam$id == t1[i],]  #Data of single tree
  c3 <- cfix[t2$src[1]]
  a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^c3)^c1, data=t2, start=list(c0=15, c1=0.6, c2=-10), control=nls.lm.control(maxiter=999), trace=F)
  a <- a0
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
#  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
}
res[["LnFixC3"]] <- resCoef
summary(resCoef)

#h = c0*log(1 + exp(c2)*t^c3fix)^c1fix
fun <- function(data, par) {
  c1 <- par[1]
  c3 <- par[2]
  t1 <- unique(data$id)
  resSd <- rep(NA, length(t1))
  for(i in 1:length(t1)) {
    t2 <- data[data$id == t1[i],]  #Data of single tree
    a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^c3)^c1, data=t2, start=list(c0=15, c2=-10), control=nls.lm.control(maxiter=999), trace=F)
    resSd[i] <- summary(a0)$sigma
  }
  median(resSd)
}
cfix <- c(optim(par=c(0.5, 3), fun, data=datStam[datStam$src=="gut",])$par, optim(par=c(0.5, 3), fun, data=datStam[datStam$src=="jrm",])$par)
cfix <- matrix(cfix, ncol=2)
colnames(cfix) <- c("gut", "jrm")
resCoef <- resCoefTempl
t1 <- unique(datStam$id)
for(i in 1:length(t1)) {
  t2 <- datStam[datStam$id == t1[i],]  #Data of single tree
  c1 <- cfix[1,t2$src[1]]
  c3 <- cfix[2,t2$src[1]]
  a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^c3)^c1, data=t2, start=list(c0=15, c2=-10), control=nls.lm.control(maxiter=999), trace=F)
  a <- a0
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
#  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
}
res[["LnFixC1C3"]] <- resCoef
summary(resCoef)

#h = c0*log(1 + exp(c2)*t^(c2*c3fix))^c1fix
fun <- function(data, par) {
  c1 <- par[1]
  c3 <- par[2]
  t1 <- unique(data$id)
  resSd <- rep(NA, length(t1))
  for(i in 1:length(t1)) {
    t2 <- data[data$id == t1[i],]  #Data of single tree
    a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^(c2*c3))^c1, data=t2, start=list(c0=15, c2=-10), control=nls.lm.control(maxiter=999), trace=F)
    resSd[i] <- summary(a0)$sigma
  }
  median(resSd)
}
cfix <- c(optim(par=c(0.5, -0.25), fun, data=datStam[datStam$src=="gut",])$par, optim(par=c(0.5, -0.25), fun, data=datStam[datStam$src=="jrm",])$par)
cfix <- matrix(cfix, ncol=2)
colnames(cfix) <- c("gut", "jrm")
resCoef <- resCoefTempl
t1 <- unique(datStam$id)
for(i in 1:length(t1)) {
  t2 <- datStam[datStam$id == t1[i],]  #Data of single tree
  c1 <- cfix[1,t2$src[1]]
  c3 <- cfix[2,t2$src[1]]
  a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^(c2*c3))^c1, data=t2, start=list(c0=15, c2=-10), control=nls.lm.control(maxiter=999), trace=F)
  a <- a0
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
#  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
}
res[["LnFixC1C3Fun1"]] <- resCoef
summary(resCoef)

#h = c0*log(1 + exp(c2)*t^(c30fix + c2*c31fix))^c1fix
fun <- function(data, par) {
  c1 <- par[1]
  c30 <- par[2]
  c31 <- par[3]
  t1 <- unique(data$id)
  resSd <- rep(NA, length(t1))
  for(i in 1:length(t1)) {
    t2 <- data[data$id == t1[i],]  #Data of single tree
    a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^(c30 + c2*c31))^c1, data=t2, start=list(c0=15, c2=-10), control=nls.lm.control(maxiter=999), trace=F)
    resSd[i] <- summary(a0)$sigma
  }
  median(resSd)
}
cfix <- c(optim(par=c(0.5367364, 3.2058904, 0), fun, data=datStam[datStam$src=="gut",])$par, optim(par=c(0.5, 1.5, -0.12), fun, data=datStam[datStam$src=="jrm",])$par)
cfix <- matrix(cfix, ncol=2)
colnames(cfix) <- c("gut", "jrm")
resCoef <- resCoefTempl
t1 <- unique(datStam$id)
for(i in 1:length(t1)) {
  t2 <- datStam[datStam$id == t1[i],]  #Data of single tree
  c1 <- cfix[1,t2$src[1]]
  c30 <- cfix[2,t2$src[1]]
  c31 <- cfix[3,t2$src[1]]
  a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^(c30 + c2*c31))^c1, data=t2, start=list(c0=15, c2=-10), control=nls.lm.control(maxiter=999), trace=F)
  a <- a0
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
#  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
}
res[["LnFixC1C3Fun2"]] <- resCoef
summary(resCoef)

#h = c0*log(1 + exp(c2)*t^c3fix)^(c1fix/c2)
fun <- function(data, par) {
  c1 <- par[1]
  c3 <- par[2]
  t1 <- unique(data$id)
  resSd <- rep(NA, length(t1))
  for(i in 1:length(t1)) {
    t2 <- data[data$id == t1[i],]  #Data of single tree
    a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^c3)^(c1/c2), data=t2, start=list(c0=15, c2=-10), control=nls.lm.control(maxiter=999), trace=F)
    resSd[i] <- summary(a0)$sigma
  }
  median(resSd)
}
cfix <- c(optim(par=c(-6.5, 3), fun, data=datStam[datStam$src=="gut",])$par, optim(par=c(-6.5, 3), fun, data=datStam[datStam$src=="jrm",])$par)
cfix <- matrix(cfix, ncol=2)
colnames(cfix) <- c("gut", "jrm")
resCoef <- resCoefTempl
t1 <- unique(datStam$id)
for(i in 1:length(t1)) {
  t2 <- datStam[datStam$id == t1[i],]  #Data of single tree
  c1 <- cfix[1,t2$src[1]]
  c3 <- cfix[2,t2$src[1]]
  a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^c3)^(c1/c2), data=t2, start=list(c0=15, c2=-10), control=nls.lm.control(maxiter=999), trace=F)
  a <- a0
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
#  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
}
res[["LnFixC1C3Fun3"]] <- resCoef
summary(resCoef)

#h = c0*log(1 + exp(c2)*t^(c2*c3fix))^(c1fix/c2)
fun <- function(data, par) {
  c1 <- par[1]
  c3 <- par[2]
  t1 <- unique(data$id)
  resSd <- rep(NA, length(t1))
  for(i in 1:length(t1)) {
    t2 <- data[data$id == t1[i],]  #Data of single tree
    a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^(c2*c3))^(c1/c2), data=t2, start=list(c0=6, c2=-6), control=nls.lm.control(maxiter=999), trace=F)
    resSd[i] <- summary(a0)$sigma
  }
  median(resSd)
}
cfix <- c(optim(par=c(-7, -0.35), fun, lower=c(-9,-0.5),upper=c(-5,-0.2), data=datStam[datStam$src=="gut",], method="L-BFGS-B")$par, optim(par=c(-7, -0.35), fun, lower=c(-9,-0.5),upper=c(-5,-0.2), data=datStam[datStam$src=="jrm",], method="L-BFGS-B")$par)
cfix <- matrix(cfix, ncol=2)
colnames(cfix) <- c("gut", "jrm")
resCoef <- resCoefTempl
t1 <- unique(datStam$id)
for(i in 1:length(t1)) {
  t2 <- datStam[datStam$id == t1[i],]  #Data of single tree
  c1 <- cfix[1,t2$src[1]]
  c3 <- cfix[2,t2$src[1]]
  a0 <- nlsLM(h ~ c0*log(1 + exp(c2)*t^(c2*c3))^(c1/c2), data=t2, start=list(c0=15, c2=-10), control=nls.lm.control(maxiter=999), trace=F)
  a <- a0
  t3 <- with(t2, splinefun(t, h, method="monoH.FC"))
  resCoef$src[i] <- as.character(t2$src[1])
  resCoef$id[i] <- as.character(t2$id[1])
  resCoef$si[i] <- t3(100)
  resCoef$sd[i] <- summary(a)$sigma
#  resCoef[i,c("c0","c1","c2","c3")] <- coef(a)
  for(j in 1:length(ageSteps)) {
    if(max(t2$t) >= ageSteps[j]) { #nicht extrapolieren
      resCoef[i,paste("h", ageSteps[j], sep="_")] <- t3(ageSteps[j])
    }
    resCoef[i,paste("hEst", ageSteps[j], sep="_")] <- predict(a, newdata=data.frame(t=ageSteps[j]))
  }
}
res[["LnFixC1C3Fun4"]] <- resCoef

res[["LnFixC1C3Fn"]] <- res[["LnFixC1C3Fun2"]]
res[["LnFixC1C3Fun1"]] <- NULL
res[["LnFixC1C3Fun2"]] <- NULL
res[["LnFixC1C3Fun3"]] <- NULL
res[["LnFixC1C3Fun4"]] <- NULL

save(res, file="results.RData")

