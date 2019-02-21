############# Version with simplex  ########

library(rethinking)
library(xlsx)
d <- read.xlsx("rawdata.xlsx",1)
str(d)
d$group <-   as.integer(ifelse(d$test==0 & d$infected==0,1,
                        ifelse(d$test==0 & d$infected==1,2,
                        ifelse(d$test==1 & d$infected==0,3,4))))

# matrix of locations (A=1,B=2,C=3) 
Y <- array(NA,dim=c(nrow(d),11))
for (i in 1:(nrow(d))){
  Y[i,1]<-as.integer(3)
  for (j in 2:11) Y[i,j] <- as.integer(d[i,5+j]) 
}

d.stan <- list(N=nrow(d),
               species=coerce_index(d$species),
               group=d$group,
               Y=Y)

### m1: full model species, group, time #############################################

m1.code <- "
data{
  int N;
  int<lower=1,upper=3> species[N];
  int<lower=1,upper=4> group[N];
  int Y[N,11];
}
parameters{
  simplex[3] c1[3,4]; // first transition from location C by species, group 
  simplex[3] theta[3,4,9,3]; // species, group, time, location
}
model{
  for (n in 1:N){ 
    Y[n,2] ~ categorical(c1[species[n],group[n]]);  
    for (t in 3:11){
      Y[n,t] ~ categorical(theta[species[n],group[n],t-2,Y[n,t-1]]);
    }
  }
}
generated quantities{
  vector[10*N] log_lik;
  for (n in 1:N){
    log_lik[(n-1)*10+1]=categorical_lpmf(Y[n,2]|c1[species[n],group[n]]);
    for (t in 3:11){
      log_lik[(n-1)*10+t-1]=categorical_lpmf(Y[n,t]|theta[species[n],group[n],t-2,Y[n,t-1]]);
    }
  }
}"


m1 <- stan(model_code=m1.code,data=d.stan,warmup=1000,iter=3500,chains=4,cores=4)
# 88 sec

saveRDS(m1,file="m1_final.rds")
m1 <- readRDS(file="m1_final.rds")

print(m1,digits=3,pars=c("c1","theta"),prob=c(0.055,0.945))
# too long...

stan_plot(m1,pars=c("c1","theta"))
stan_trace(m1,pars=c("c1","theta"))
stan_dens(m1,pars=c("c1","theta"))

WAIC(m1)
# 5303.5 (98.3)


### m2 simplex: drop time  #########

m2.code <- "
data{
  int N;
  int<lower=1,upper=3> species[N];
  int<lower=1,upper=4> group[N];
  int Y[N,11];
}
parameters{
  simplex[3] c1[3,4]; // species, group
  simplex[3] theta[3,4,3]; // species, group, category
}
model{
  for (n in 1:N){ 
    Y[n,2] ~ categorical(c1[species[n],group[n]]);  
    for (t in 3:11){
      Y[n,t] ~ categorical(theta[species[n],group[n],Y[n,t-1]]);
    }
  }
}
generated quantities{
  vector[10*N] log_lik;
  for (n in 1:N){
    log_lik[(n-1)*10+1]=categorical_lpmf(Y[n,2]|c1[species[n],group[n]]);
    for (t in 3:11){
      log_lik[(n-1)*10+t-1]=categorical_lpmf(Y[n,t]|theta[species[n],group[n],Y[n,t-1]]);
    }
  }
}"


m2 <- stan(model_code=m2.code,data=d.stan,warmup=1000,iter=3500,chains=4,cores=4)
# 71 sec

saveRDS(m2,file="m2_final.rds")
m2 <- readRDS(file="m2_final.rds")

print(m2,digits=3,pars=c("c1","theta"),prob=c(0.055,0.945))

stan_plot(m2,pars=c("c1","theta"))
stan_trace(m2,pars=c("c1","theta"))
stan_dens(m2,pars=c("c1","theta"))

WAIC(m2)
# m1: 5303.6 (98.3)
# m2: 4995.6 (109.8)
# dropping time improves model

### /m2 simplex: drop time  #########



### m3 simplex: drop time and group  #########

m3.code <- "
data{
  int N;
  int<lower=1,upper=3> species[N];
  int<lower=1,upper=4> group[N];
  int Y[N,11];
}
parameters{
  simplex[3] c1[3]; // species
  simplex[3] theta[3,3]; // species, category
}
model{
  for (n in 1:N){ 
    Y[n,2] ~ categorical(c1[species[n]]);  
    for (t in 3:11){
      Y[n,t] ~ categorical(theta[species[n],Y[n,t-1]]);
    }
  }
}
generated quantities{
  vector[10*N] log_lik;
  for (n in 1:N){
    log_lik[(n-1)*10+1]=categorical_lpmf(Y[n,2]|c1[species[n]]);
    for (t in 3:11){
      log_lik[(n-1)*10+t-1]=categorical_lpmf(Y[n,t]|theta[species[n],Y[n,t-1]]);
    }
  }
}"


m3 <- stan(model_code=m3.code,data=d.stan,warmup=1000,iter=3500,chains=4,cores=4)
# 17 sec

saveRDS(m3,file="m3_final.rds")
# m3 <- readRDS(file="m3_final.rds")

print(m3,digits=3,pars=c("c1","theta"),prob=c(0.055,0.945))

stan_plot(m3,pars=c("c1","theta"))
stan_trace(m3,pars=c("c1","theta"))
stan_dens(m3,pars=c("c1","theta"))

WAIC(m3)
# m1: 5303.6 (98.3)
# m2: 4995.6 (109.8)
# dropping time improves model wrt m1
# m3: 5138.0 (112.4)
# dropping group does not improve wrt m2

### /m3 simplex: drop time and group  #########



### m4 simplex: drop species  #########

m4.code <- "
data{
  int N;
  int<lower=1,upper=3> species[N];
  int<lower=1,upper=4> group[N];
  int Y[N,11];
}
parameters{
  simplex[3] c1[4]; // group
  simplex[3] theta[4,9,3]; // group, time, category
}
model{
  for (n in 1:N){ 
    Y[n,2] ~ categorical(c1[group[n]]);  
    for (t in 3:11){
      Y[n,t] ~ categorical(theta[group[n],t-2,Y[n,t-1]]);
    }
  }
}
generated quantities{
  vector[10*N] log_lik;
  for (n in 1:N){
    log_lik[(n-1)*10+1]=categorical_lpmf(Y[n,2]|c1[group[n]]);
    for (t in 3:11){
      log_lik[(n-1)*10+t-1]=categorical_lpmf(Y[n,t]|theta[group[n],t-2,Y[n,t-1]]);
    }
  }
}"


m4 <- stan(model_code=m4.code,data=d.stan,warmup=1000,iter=3500,chains=4,cores=4)
# 42 sec

saveRDS(m4,file="m4_final.rds")
# m4 <- readRDS(file="m4_final.rds")

print(m4,digits=3,pars=c("c1","theta"),prob=c(0.055,0.945))

stan_plot(m4,pars=c("c1","theta"))
stan_trace(m4,pars=c("c1","theta"))
stan_dens(m4,pars=c("c1","theta"))

WAIC(m4)
# m1: 5303.6 (98.3)
# m2: 4995.6 (109.8)
# dropping time improves model wrt m1
# m3: 5138.0 (112.4)
# dropping group does not improve wrt m2
# m4: 5244.5 (109.5)
# dropping species improves wrt m1

### /m4 simplex: drop species  #########


### m5 simplex: drop time and species  #########

m5.code <- "
data{
  int N;
  int<lower=1,upper=3> species[N];
  int<lower=1,upper=4> group[N];
  int Y[N,11];
}
parameters{
  simplex[3] c1[4]; // group
  simplex[3] theta[4,3]; // group, category
}
model{
  for (n in 1:N){ 
    Y[n,2] ~ categorical(c1[group[n]]);  
    for (t in 3:11){
      Y[n,t] ~ categorical(theta[group[n],Y[n,t-1]]);
    }
  }
}
generated quantities{
  vector[10*N] log_lik;
  for (n in 1:N){
    log_lik[(n-1)*10+1]=categorical_lpmf(Y[n,2]|c1[group[n]]);
    for (t in 3:11){
      log_lik[(n-1)*10+t-1]=categorical_lpmf(Y[n,t]|theta[group[n],Y[n,t-1]]);
    }
  }
}"


m5 <- stan(model_code=m5.code,data=d.stan,warmup=1000,iter=3500,chains=4,cores=4)


saveRDS(m5,file="m5_final.rds")
# m5 <- readRDS(file="m5_final.rds")

print(m5,digits=3,pars=c("c1","theta"),prob=c(0.055,0.945))

stan_plot(m5,pars=c("c1","theta"))
stan_trace(m5,pars=c("c1","theta"))
stan_dens(m5,pars=c("c1","theta"))

WAIC(m5)
# m1: 5303.6 (98.3)
# m2: 4995.6 (109.8)
# dropping time improves model wrt m1
# m3: 5138.0 (112.4)
# dropping group does not improve wrt m2
# m4: 5244.5 (109.5)
# dropping species improves wrt m1
# m5: 5100.6 (111.7)
# dropping time improves wrt m4
# but dropping species does not improve wrt m2


### /m5 simplex: drop time and species  #########


### m6 simplex: drop group  #########

m6.code <- "
data{
  int N;
  int<lower=1,upper=3> species[N];
  int<lower=1,upper=4> group[N];
  int Y[N,11];
}
parameters{
  simplex[3] c1[3]; // species, group
  simplex[3] theta[3,9,3]; // species, group, time, category
}
model{
  for (n in 1:N){ 
    Y[n,2] ~ categorical(c1[species[n]]);  
    for (t in 3:11){
      Y[n,t] ~ categorical(theta[species[n],t-2,Y[n,t-1]]);
    }
  }
}
generated quantities{
  vector[10*N] log_lik;
  for (n in 1:N){
    log_lik[(n-1)*10+1]=categorical_lpmf(Y[n,2]|c1[species[n]]);
    for (t in 3:11){
      log_lik[(n-1)*10+t-1]=categorical_lpmf(Y[n,t]|theta[species[n],t-2,Y[n,t-1]]);
    }
  }
}"


m6 <- stan(model_code=m6.code,data=d.stan,warmup=1000,iter=3500,chains=4,cores=4)

saveRDS(m6,file="m6_final.rds")
# m6 <- readRDS(file="m6_final.rds")

print(m6,digits=3,pars=c("c1","theta"),prob=c(0.055,0.945))

stan_plot(m6,pars=c("c1","theta"))
stan_trace(m6,pars=c("c1","theta"))
stan_dens(m6,pars=c("c1","theta"))

WAIC(m6)
# m1: 5303.6 (98.3)
# m2: 4995.6 (109.8)
# dropping time improves model wrt m1
# m3: 5138.0 (112.4)
# dropping group does not improve wrt m2
# m4: 5244.5 (109.5)
# dropping species improves wrt m1
# m5: 5100.6 (111.7) full minus species & time (group left)
# improves wrt m4 but not wrt m2
# m6: 5269.2 (111.3) full minus group (species & time left)
# improves wrt full but not wrt m3

### /m6 simplex: drop group  #########


### m7 simplex: drop group & species #########

m7.code <- "
data{
  int N;
  int<lower=1,upper=3> species[N];
  int<lower=1,upper=4> group[N];
  int Y[N,11];
}
parameters{
  simplex[3] c1; // 
  simplex[3] theta[9,3]; //time, category
}
model{
  for (n in 1:N){ 
    Y[n,2] ~ categorical(c1);  
    for (t in 3:11){
      Y[n,t] ~ categorical(theta[t-2,Y[n,t-1]]);
    }
  }
}
generated quantities{
  vector[10*N] log_lik;
  for (n in 1:N){
    log_lik[(n-1)*10+1]=categorical_lpmf(Y[n,2]|c1);
    for (t in 3:11){
      log_lik[(n-1)*10+t-1]=categorical_lpmf(Y[n,t]|theta[t-2,Y[n,t-1]]);
    }
  }
}"


m7 <- stan(model_code=m7.code,data=d.stan,warmup=1000,iter=3500,chains=4,cores=4)
# 20 sec [long plus some bad transitions]

saveRDS(m7,file="m7_final.rds")
# m7 <- readRDS(file="m7_final.rds")

print(m7,digits=3,pars=c("c1","theta"),prob=c(0.055,0.945))

stan_plot(m7,pars=c("c1","theta"))
stan_trace(m7,pars=c("c1","theta"))
stan_dens(m7,pars=c("c1","theta"))

WAIC(m7)
# m1: 5303.6 (98.3)
# m2: 4995.6 (109.8)
# dropping time improves model wrt m1
# m3: 5138.0 (112.4)
# dropping group does not improve wrt m2
# m4: 5244.5 (109.5)
# dropping species improves wrt m1
# m5: 5100.6 (111.7) full minus species & time (group left)
# improves wrt m4 but not wrt m2
# m6: 5269.2 (111.3) full minus group (species & time left)
# improves wrt full but not wrt m3
# m7: 5288.5 (113.8) full minus group & species (only time left)
# improves wrt m1 but not m4, m6
# Conclusion: m2 is best, just drop time

### /m7 simplex: drop group & species #########

p2 <- extract(m2,permute=2)
str(p2)

nsim <- 10000
nsp <- 3
ngr <- 4

# Posterior (species, group) for prob to be at time (2..11) and loc (1..3)
P <- array(dim=c(nsim,nsp,ngr,10,3))

for (sp in 1:nsp) for (gr in 1:ngr) for (loc in 1:3){
  P[,sp,gr,1,loc] <- p2$c1[,sp,gr,loc]
}

for (sp in 1:nsp) for (gr in 1:ngr) for (time in 2:10) for (loc in 1:3){
  P[,sp,gr,time,loc] <- p2$theta[,sp,gr,1,loc]*P[,sp,gr,time-1,1]+
                        p2$theta[,sp,gr,2,loc]*P[,sp,gr,time-1,2]+
                        p2$theta[,sp,gr,3,loc]*P[,sp,gr,time-1,3]
}

# Posterior probs to be in locs (by spec & group)

Q <- array(0,dim=c(nsim,nsp,ngr,3))
for (sp in 1:nsp) for (gr in 1:ngr) for (loc in 1:3){ 
  for (time in 1:10){
    Q[,sp,gr,loc] <- Q[,sp,gr,loc]+P[,sp,gr,time,loc]
  }
  Q[,sp,gr,loc] <- Q[,sp,gr,loc]/10
}

# Data for coefficient plots
PAminB.mean <- apply(Q[,,,1]-Q[,,,2],c(2,3),mean)
PAminB.min  <- apply(Q[,,,1]-Q[,,,2],c(2,3),HPDI)[1,,]
PAminB.max  <- apply(Q[,,,1]-Q[,,,2],c(2,3),HPDI)[2,,]

PAminC.mean <- apply(Q[,,,1]-Q[,,,3],c(2,3),mean)
PAminC.min  <- apply(Q[,,,1]-Q[,,,3],c(2,3),HPDI)[1,,]
PAminC.max  <- apply(Q[,,,1]-Q[,,,3],c(2,3),HPDI)[2,,]

PBminC.mean <- apply(Q[,,,2]-Q[,,,3],c(2,3),mean)
PBminC.min  <- apply(Q[,,,2]-Q[,,,3],c(2,3),HPDI)[1,,]
PBminC.max  <- apply(Q[,,,2]-Q[,,,3],c(2,3),HPDI)[2,,]

PAminB.mean
as.vector(PAminB.mean) # stacked columns (group)

ylabs <- c("UC","IC","UT","IT")
species <- c("Species EB","Species GF","Species GP")
dgAminB <- data.frame(x=as.vector(t(PAminB.mean)),
                      xmin=as.vector(t(PAminB.min)),
                      xmax=as.vector(t(PAminB.max)),
                      y=rep(1:4,3),
                      Species=as.factor(rep(species,each=4)))

library(ggplot2)
ggplot(data=dgAminB,aes(x=x,y=y))+
geom_vline(xintercept=0,col="grey")+
geom_point(size=3)+
facet_grid(Species~.)+
theme_bw()+
geom_errorbarh(aes(xmin=xmin,xmax=xmax),height=0,lwd=1)+
scale_x_continuous(limits=c(-0.2,0.3),breaks=seq(-0.2,0.3,0.1))+
scale_y_continuous(limits=c(0.5,4.5),breaks=1:4,labels=ylabs)+
labs(x="\nP(A) - P(B)",y="Category\n")



dgAminC <- data.frame(x=as.vector(t(PAminC.mean)),
                      xmin=as.vector(t(PAminC.min)),
                      xmax=as.vector(t(PAminC.max)),
                      y=rep(1:4,3),
                      Species=as.factor(rep(species,each=4)))

library(ggplot2)
ggplot(data=dgAminC,aes(x=-x,y=y))+
geom_point(size=3)+
facet_grid(Species~.)+
theme_bw()+
geom_errorbarh(aes(xmin=-xmin,xmax=-xmax),height=0,lwd=1)+
scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
scale_y_continuous(limits=c(0.5,4.5),breaks=1:4,labels=ylabs)+
labs(x="\nP(C) - P(A)",y="Category\n")



dgCminB <- data.frame(x=-as.vector(t(PBminC.mean)),
                      xmin=-as.vector(t(PBminC.min)),
                      xmax=-as.vector(t(PBminC.max)),
                      y=rep(1:4,3),
                      Species=as.factor(rep(species,each=4)))

library(ggplot2)
ggplot(data=dgCminB,aes(x=x,y=y))+
geom_point(size=3)+
facet_grid(Species~.)+
theme_bw()+
geom_errorbarh(aes(xmin=xmin,xmax=xmax),height=0,lwd=1)+
scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
scale_y_continuous(limits=c(0.5,4.5),breaks=1:4,labels=ylabs)+
labs(x="\nP(C) - P(B)",y="Category\n")

# /coefficient plots


## new time plots

# data
PAT.mean <- apply(P[,,,,1],c(2,3,4),mean)
as.vector(PAT.mean)
PAT.min <- apply(P[,,,,1],c(2,3,4),HPDI)[1,,,]
PAT.max <- apply(P[,,,,1],c(2,3,4),HPDI)[2,,,]

PBT.mean <- apply(P[,,,,2],c(2,3,4),mean)
PBT.min <- apply(P[,,,,2],c(2,3,4),HPDI)[1,,,]
PBT.max <- apply(P[,,,,2],c(2,3,4),HPDI)[2,,,]

PCT.mean <- apply(P[,,,,3],c(2,3,4),mean)
PCT.min <- apply(P[,,,,3],c(2,3,4),HPDI)[1,,,]
PCT.max <- apply(P[,,,,3],c(2,3,4),HPDI)[2,,,]

species <- c("Species EB","Species GF","Species GP")
groups <- c("UC","IC","UT","IT")
dgT <- data.frame(time=30*rep(1:10,each=12),
                  Species=rep(species,40),
                  Group=rep(groups,each=3,times=10),
                  yA=as.vector(PAT.mean),
                  yAmin=as.vector(PAT.min),
                  yAmax=as.vector(PAT.max),
                  yB=as.vector(PBT.mean),
                  yBmin=as.vector(PBT.min),
                  yBmax=as.vector(PBT.max),
                  yC=as.vector(PCT.mean),
                  yCmin=as.vector(PCT.min),
                  yCmax=as.vector(PCT.max))

library(ggplot2)
ggplot(data=dgT,aes(x=time,y=yA))+
geom_line(col="red")+
geom_line(aes(x=time,y=yB),col="blue")+
geom_line(aes(x=time,y=yC),col="green")+
geom_ribbon(aes(ymin=yAmin,ymax=yAmax),alpha=0.2,fill="red")+
geom_ribbon(aes(ymin=yBmin,ymax=yBmax),alpha=0.2,fill="blue")+
geom_ribbon(aes(ymin=yCmin,ymax=yCmax),alpha=0.2,fill="green")+
facet_grid(Species~Group)+
theme_bw()+
scale_x_continuous(limits=c(0,300),breaks=seq(0,300,100))+
scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
labs(x="\nTime (s)",y="Proportion present\n")

## raw data time plots for comparison

ns <- table(d$species,d$group)

nsp <- 3
ngr <- 4

PATraw <- array(dim=c(nsp,ngr,10))
PBTraw <- array(dim=c(nsp,ngr,10))
PCTraw <- array(dim=c(nsp,ngr,10))

for (i in 1:10) PATraw[,,i]<-tapply(d[,6+i],list(d$species,d$group),function(x) sum(x=="A"))/ns
for (i in 1:10) PBTraw[,,i]<-tapply(d[,6+i],list(d$species,d$group),function(x) sum(x=="B"))/ns
for (i in 1:10) PCTraw[,,i]<-tapply(d[,6+i],list(d$species,d$group),function(x) sum(x=="C"))/ns

# heh heh, nice code...

dgTraw <- data.frame(time=30*rep(1:10,each=12),
                    Species=rep(species,40),
                    Group=rep(groups,each=3,times=10),
                    yA=as.vector(PATraw),
                    yB=as.vector(PBTraw),
                    yC=as.vector(PCTraw))

ggplot(data=dgTraw,aes(x=time,y=yA))+
geom_line(col="red")+
geom_line(aes(x=time,y=yB),col="blue")+
geom_line(aes(x=time,y=yC),col="green")+
facet_grid(Species~Group)+
theme_bw()+
scale_x_continuous(limits=c(0,300),breaks=seq(0,300,100))+
scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
labs(x="\nTime (s)",y="Proportion present\n")


# combine raw with smoothed

species <- c("Species EB","Species GF","Species GP")
groups <- c("UC","IC","UT","IT")
dgT <- data.frame(time=30*rep(1:10,each=12),
                  Species=rep(species,40),
                  Group=rep(groups,each=3,times=10),
                  yA=as.vector(PAT.mean),
                  yAmin=as.vector(PAT.min),
                  yAmax=as.vector(PAT.max),
                  yB=as.vector(PBT.mean),
                  yBmin=as.vector(PBT.min),
                  yBmax=as.vector(PBT.max),
                  yC=as.vector(PCT.mean),
                  yCmin=as.vector(PCT.min),
                  yCmax=as.vector(PCT.max),
                  yAraw=as.vector(PATraw),
                  yBraw=as.vector(PBTraw),
                  yCraw=as.vector(PCTraw))

library(ggplot2)
ggplot(data=dgT,aes(x=time,y=yA))+
geom_line(col="darkred",lwd=1)+
geom_line(aes(x=time,y=yB),col="darkblue",lwd=1)+
geom_line(aes(x=time,y=yC),col="darkgreen",lwd=1)+
geom_line(aes(x=time,y=yAraw),col="red")+
geom_line(aes(x=time,y=yBraw),col="blue")+
geom_line(aes(x=time,y=yCraw),col="green")+
geom_ribbon(aes(ymin=yAmin,ymax=yAmax),alpha=0.2,fill="red")+
geom_ribbon(aes(ymin=yBmin,ymax=yBmax),alpha=0.2,fill="blue")+
geom_ribbon(aes(ymin=yCmin,ymax=yCmax),alpha=0.2,fill="green")+
facet_grid(Species~Group)+
theme_bw()+
scale_x_continuous(limits=c(0,300),breaks=seq(0,300,100))+
scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
labs(x="\nTime (s)",y="Proportion present\n")

# That fits very nicely!!!


### Model comparison


weight <- function(x){
  d <- x-min(x)
  log_denom <- log_sum_exp(-0.5*d)
  logw <- -0.5*d-log_denom
  return(c(d,exp(logw)))
}

x <- c(5303.6,4995.6,5138.0,5244.5,5100.6,5269.2,5288.5)
print(round(weight(x),digits=2))

