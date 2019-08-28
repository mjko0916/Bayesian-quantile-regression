library(rjags)

set.seed(200)

##### bus data #####
bus <- read.csv("final.csv")[,-1]
bus <- bus[sample(1:nrow(bus), 2000),]
x = as.matrix(bus[,c(6,7,8,12,13,16,17,20,21,22,24:30)])
for(i in 1:ncol(x)) x[,i]=(x[,i]-mean(x[,i]))/sd(x[,i])
x = cbind(rep(1,nrow(x)), x)
x = data.matrix(x)
colnames(x)[1] = "intercept"
head(x)
y = bus[,10]
K = ncol(x)


################## JAGS_BQR_GVS #################
library(quantreg)
require(plyr); require(data.table)
codaSamples <- list()
beta.samples <- list(); gamma.samples <- list()
nAdapt = 10000; nUpdate = 10000; nIter = 20000; nChains = 3; n.thin=10
p.set = seq(0.1, 0.9, length = 9)



for(ip in 1:9) {
  p=p.set[ip]
  rq.out=rq.fit.br(x,y,p, ci=T)
  pseudo.mean.beta=rq.out$coefficients[1:K,1]
  pseudo.sd.beta=(rq.out$coefficients[1:K,3]-rq.out$coefficients[1:K,2])/(2*1.96)
  pseudo.var.beta=pseudo.sd.beta^2
  z=rep(1,length(y))
  dataList=list(p=p, K=K, y=y, x=x, z=z, 
                pseudo.mean.beta=pseudo.mean.beta,
                pseudo.var.beta=pseudo.var.beta)
  gammaInit=rep(0,K)
  initsList=list(beta=pseudo.mean.beta, gamma=gammaInit) 
  
  jagsModel=jags.model(file="model_BQR_GVS.txt", data=dataList,
                       inits=initsList,
                       n.chains=nChains, n.adapt=nAdapt)
  update(jagsModel, n.iter=nUpdate)
  codaSamples[[ip]]=coda.samples(jagsModel, variable.names=c("beta","gamma"), thin=n.thin,
                                 n.iter=nIter)
  
  mcmcSamples=as.matrix(codaSamples[[ip]])
  beta.samples[[ip]]=mcmcSamples[,1:K]
  gamma.samples[[ip]]=mcmcSamples[,(K+1):(K+K)]
}

### beta 수렴 진단
conv <- c()
for (i in 2:9){ 
  conv[i] <- gelman.diag(codaSamples[[i]][,1:K])
}

# ACF plots
par(mfrow = c(2,3))
for (i in 1:length(p.set)) {
  for (j in 2:K) acf(beta.samples[[i]][,j], main = colnames(x)[j])
} 
par(mfrow=c(1,1))
for (j in 2:K) acf(beta.samples[[9]][,j], main = colnames(x)[j])

### beta 추정치
m <- list(); gamma.hat <- list(); beta.selected.hat <- list()
for (ip in 2:length(p.set)){
  m=gamma.samples[[ip]]
  mm=as.data.table(m)[, .N, by = eval(paste0("gamma[", seq_len(ncol(m)), "]"))]
  colnames(mm)=c("g0","g1","g2","g3","g4","g5","g6","g7","g8","g9","g10","g11","g12","g13","g14","g15","g16","g17","N")
  mm.order=order(mm$N, decreasing=T)
  mm$N=round( mm$N/(nIter*nChains/n.thin),4)
  gamma.hat[[ip]]=as.numeric(mm[which.max(mm$N)])
  
  gamma.samples.collapsed <- apply(m, 1, function(x) paste(x, collapse=" "))
  gamma.hat.collapsed=paste(gamma.hat[[ip]][1:K], collapse=" ")
  id.selected=which(gamma.samples.collapsed==gamma.hat.collapsed)
  
  beta.samples.selected=beta.samples[[ip]][id.selected,]
  colnames(beta.samples.selected)=c("b0","b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14","b15","b16","b17")
  beta.samples.selected2=beta.samples.selected[,gamma.hat[[ip]][1:K]==1]
  beta.selected.hat[[ip]]=apply(beta.samples.selected2,2,function(x)quantile(x, c(0.025,0.5, 0.975)))
}

### beta hat plot

beta.hat = beta.cl=beta.cu= matrix(0, 9,K)
for(ip in 1:length(p.set)) {
  mcmcSamples <- as.matrix(codaSamples[[ip]])
  
  beta.hat[(ip), 1:K] = apply( mcmcSamples[,1:K],2,quantile, 0.5)
  beta.cl[(ip), 1:K] = apply( mcmcSamples[,1:K],2,quantile, 0.025)
  beta.cu[(ip), 1:K] = apply( mcmcSamples[,1:K],2,quantile, 0.975)
}

library(RColorBrewer)

colpal <- brewer.pal(n=9, name="YlOrBr")
library(tidyverse)

beta.hat <- data.frame(beta.hat)
beta.hat$p <- p.set
colnames(beta.hat)[1:18] <- colnames(x)

beta.cl <- data.frame(beta.cl)
beta.cl$p <- p.set
colnames(beta.cl)[1:18] <- colnames(x)

beta.cu <- data.frame(beta.cu)
beta.cu$p <- p.set
colnames(beta.cu)[1:18] <- colnames(x)


gg<-list()
for(i in 2:K){
  gg1 <- eval(substitute(ggplot() + geom_line(data=beta.hat, aes(x=p.set[2:9], y=beta.hat[, i]), col=colpal[6], size=1.5) + 
                           geom_line(data=beta.cl, aes(x=p.set[2:9], y=beta.cl[, i]), col=colpal[8], size=1.5) +
                           geom_line(data=beta.cu, aes(x=p.set[2:9], y=beta.cu[, i]), col=colpal[8], size=1.5) +
                           geom_hline(yintercept=0, linetype="dashed") + theme_bw() +
                           labs(x="p", y="") + ggtitle(paste(colnames(x)[i], "-", "승차")) +
                           theme(plot.title = element_text(size=20)) +
                           scale_x_continuous(breaks=seq(0.1,0.9,0.1), labels=seq(0.1,0.9,0.1)), 
                         list(i=i)))
  
  gg[[i-1]] <- gg1
}
library(gridExtra)
grid.arrange(gg[[1]], gg[[2]], gg[[3]], gg[[4]], gg[[5]], gg[[6]], ncol=3, nrow=2)
grid.arrange(gg[[7]], gg[[8]], gg[[9]], gg[[10]], gg[[11]], gg[[12]], ncol=3, nrow=2)
grid.arrange(gg[[13]], gg[[14]], gg[[15]], gg[[16]], gg[[17]], ncol=3, nrow=2)
