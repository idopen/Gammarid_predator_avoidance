library(tidyverse)
library(readxl)
library(cowplot)

d <- read_excel("rawdata_2.xlsx",1)
str(d)

d <- as_tibble(d)

### Turn into "long" format:
d2 <- d %>% gather(`2`:`11`,key = "time", value = "location")
str(d2)

## Time from character to integer and then sort by individual and time:
d2$time <- as.integer(d2$time)
d2 <- d2 %>% arrange(individual,time)

## Locations to integers (A=1,B=2,C=3) and call it y
d2$location <- ifelse(d2$location == "A",1,ifelse(d2$location == "B",2,3))
d2$location <- as.integer(d2$location)
d2 <- rename(d2, y = location)

## Create variable for location at previous time
## Except for time=2 because we know at time=1 y=3
d2$y_prev <- 3
for (i in 1:nrow(d2)){
  if (d2$time[i] > 2) d2$y_prev[i] <- d2$y[i-1]
}

## Make factors
d2$y_prev <- as.factor(d2$y_prev)
d2$time <- as.factor(d2$time)
d2$species <- as.factor(d2$species)
d2$test <- as.factor(d2$test)
d2$individual <- as.factor(d2$individual)
d2$infected <- as.factor(d2$infected)


####### brms models ###################################################################

library(brms)

d2$y1 <- d2$y
d2$y2 <- d2$y
d2$pred1 <- d2$species:d2$infected:d2$test
d2$pred2 <- d2$pred1:d2$y_prev:d2$time
d2$sub1 <- ifelse(d2$time=="2",1,0)
d2$sub2 <- ifelse(d2$sub1==0,1,0)

bform1 <- bf(y1 | subset(sub1) ~ pred1) + categorical() +
  bf(y2 | subset(sub2) ~ pred2) + categorical() + set_rescor(F)

m1 <- brm(bform1, 
          prior = prior(normal(0, 5), class = "b"),
          data = d2,
          warmup = 1000,
          iter = 3500,
          chains = 4,
          cores = 4,
          file = "m1_brms")

## 3900-4000 sec

summary(m1)


## Drop time from model

d2$pred2 <- d2$pred1:d2$y_prev

bform2 <- bf(y1 | subset(sub1) ~ pred1) + categorical() +
  bf(y2 | subset(sub2) ~ pred2) + categorical() + set_rescor(F)

m2 <- brm(bform2, 
          prior = prior(normal(0, 5), class = "b"),
          data = d2,
          warmup = 1000,
          iter = 3500,
          chains = 4,
          cores = 4,
          file = "m2_brms")

## 490-530 sec

## diagnostics:
launch_shinystan(m2)
## all looks good

## Drop time and infected:test

d2$pred1 <- d2$species
d2$pred2 <- d2$pred1:d2$y_prev

bform3 <- bf(y1 | subset(sub1) ~ pred1) + categorical() +
  bf(y2 | subset(sub2) ~ pred2) + categorical() + set_rescor(F)

m3 <- brm(bform3, 
          prior = prior(normal(0, 5), class = "b"),
          data = d2,
          warmup = 1000,
          iter = 3500,
          chains = 4,
          cores = 4,
          file = "m3_brms")

## 180 sec

## m4: drop species from m1

d2$pred1 <- d2$infected:d2$test
d2$pred2 <- d2$pred1:d2$y_prev:d2$time

bform4 <- bf(y1 | subset(sub1) ~ pred1) + categorical() +
  bf(y2 | subset(sub2) ~ pred2) + categorical() + set_rescor(F)

m4 <- brm(bform4, 
          prior = prior(normal(0, 5), class = "b"),
          data = d2,
          warmup = 1000,
          iter = 3500,
          chains = 4,
          cores = 4,
          file = "m4_brms")

## 1050 sec

## m5: drop time and species

d2$pred1 <- d2$infected:d2$test
d2$pred2 <- d2$pred1:d2$y_prev

bform5 <- bf(y1 | subset(sub1) ~ pred1) + categorical() +
  bf(y2 | subset(sub2) ~ pred2) + categorical() + set_rescor(F)

m5 <- brm(bform5, 
          prior = prior(normal(0, 5), class = "b"),
          data = d2,
          warmup = 1000,
          iter = 3500,
          chains = 4,
          cores = 4,
          file = "m5_brms")

## 260 sec

## m6: drop infected:test

d2$pred1 <- d2$species
d2$pred2 <- d2$pred1:d2$y_prev:d2$time

bform6 <- bf(y1 | subset(sub1) ~ pred1) + categorical() +
  bf(y2 | subset(sub2) ~ pred2) + categorical() + set_rescor(F)

m6 <- brm(bform6, 
          prior = prior(normal(0, 5), class = "b"),
          data = d2,
          warmup = 1000,
          iter = 3500,
          chains = 4,
          cores = 4,
          file = "m6_brms")

## 820-840 sec

## m7: drop infected:test:species

d2$pred2 <- d2$y_prev:d2$time

bform7 <- bf(y1 | subset(sub1) ~ 1) + categorical() +
  bf(y2 | subset(sub2) ~ pred2) + categorical() + set_rescor(F)

m7 <- brm(bform7, 
          prior = prior(normal(0, 5), class = "b"),
          data = d2,
          warmup = 1000,
          iter = 3500,
          chains = 4,
          cores = 4,
          file = "m7_brms")

## 210-260 sec

######### Model comparison ###################

waics <- array(dim=c(7,2)) # 
waics[1,1] <- WAIC(m1,resp = "y1")$waic + WAIC(m1, resp = "y2")$waic
waics[1,2] <- sqrt(WAIC(m1,resp = "y1")$se_waic^2 + 
                     WAIC(m1, resp = "y2")$se_waic^2)
waics[2,1] <- WAIC(m2,resp = "y1")$waic + WAIC(m2, resp = "y2")$waic
waics[2,2] <- sqrt(WAIC(m2,resp = "y1")$se_waic^2 + 
                     WAIC(m2, resp = "y2")$se_waic^2)
waics[3,1] <- WAIC(m3,resp = "y1")$waic + WAIC(m3, resp = "y2")$waic
waics[3,2] <- sqrt(WAIC(m3,resp = "y1")$se_waic^2 + 
                     WAIC(m3, resp = "y2")$se_waic^2)
waics[4,1] <- WAIC(m4,resp = "y1")$waic + WAIC(m4, resp = "y2")$waic
waics[4,2] <- sqrt(WAIC(m4,resp = "y1")$se_waic^2 + 
                     WAIC(m4, resp = "y2")$se_waic^2)
waics[5,1] <- WAIC(m5,resp = "y1")$waic + WAIC(m5, resp = "y2")$waic
waics[5,2] <- sqrt(WAIC(m5,resp = "y1")$se_waic^2 + 
                     WAIC(m5, resp = "y2")$se_waic^2)
waics[6,1] <- WAIC(m6,resp = "y1")$waic + WAIC(m6, resp = "y2")$waic
waics[6,2] <- sqrt(WAIC(m6,resp = "y1")$se_waic^2 + 
                     WAIC(m6, resp = "y2")$se_waic^2)
waics[7,1] <- WAIC(m7,resp = "y1")$waic + WAIC(m7, resp = "y2")$waic
waics[7,2] <- sqrt(WAIC(m7,resp = "y1")$se_waic^2 + 
                     WAIC(m7, resp = "y2")$se_waic^2)
saveRDS(waics,"waics.rds")
waics <- readRDS("waics.rds")

waic_weights <- function(x){
  delta <- x-min(x)
  log_denom <- rethinking::log_sum_exp(-0.5*delta)
  logw <- -0.5*delta - log_denom
  return(c(delta,exp(logw)))
}

round(waic_weights(waics[,1]),2)

## Model m2 wins.

## We continue with model m2 and extract the posterior samples

d2$pred1 <- d2$species:d2$infected:d2$test
d2$pred2 <- d2$pred1:d2$y_prev

d.pred1 <- expand.grid(pred1 = levels(d2$pred1), sub1 = 1) 
d.pred2 <- expand.grid(pred2 = levels(d2$pred2), sub2 = 1)

## Posterior samples:
fit2_y1 <- fitted(m2, resp = "y1", newdata = d.pred1, summary = F)
fit2_y2 <- fitted(m2, resp = "y2", newdata = d.pred2, summary = F)

## Store samples in a more ordered way in array:
n_sim <- 10000
n_species <- 3
n_treat <- 4
n_time <- 10
n_loc <- 3

## Array containing samples of posterior probabilities to be in a 
## specific location at a specific time, depending on species and treatment
## Calculated from formula [3] in the manuscript
P <- array(dim=c(n_sim,n_species,n_treat,n_time,n_loc))

## species: EB, GF, GP
## treat (infected, test) 0:0,0:1,1:0,1:1

## For time=1, submodel 1
for (species in 1:n_species) for (treat in 1:n_treat) for (loc in 1:n_loc){
  index <- (species-1)*n_treat + treat
  P[,species,treat,1,loc] <- fit2_y1[,index,loc]
}

## For time>1, submodel 2, using formula [3] from manuscript
for (species in 1:n_species) for (treat in 1:n_treat) 
  for (time in 2:10) for (loc in 1:n_loc){
  index <- (species-1)*n_treat*n_loc + (treat-1)*n_loc
  P[,species,treat,time,loc] <- fit2_y2[,index+1,loc]*P[,species,treat,time-1,1]+
    fit2_y2[,index+2,loc]*P[,species,treat,time-1,2]+
    fit2_y2[,index+3,loc]*P[,species,treat,time-1,3]
}

## Array containing samples of posterior probabilities
## to be in a given location (averaged over time),
## calculated according to formula [4] in manuscript

Q <- array(0,dim=c(n_sim,n_species,n_treat,n_loc))

for (species in 1:n_species) for (treat in 1:n_treat)
  for (loc in 1:n_loc) for (time in 1:n_time){
    Q[,species,treat,loc] <- Q[,species,treat,loc] + P[,species,treat,time,loc]/n_time
  }


## Figure 2 #########################################################

## Calculate time-specific model-fitted probabilities to be in location A
## mean and HPDI:
PAT.mean <- apply(P[,,,,1],c(2,3,4),mean)
PAT.min <- apply(P[,,,,1],c(2,3,4),rethinking::HPDI)[1,,,]
PAT.max <- apply(P[,,,,1],c(2,3,4),rethinking::HPDI)[2,,,]

## Ditto location B
PBT.mean <- apply(P[,,,,2],c(2,3,4),mean)
PBT.min <- apply(P[,,,,2],c(2,3,4),rethinking::HPDI)[1,,,]
PBT.max <- apply(P[,,,,2],c(2,3,4),rethinking::HPDI)[2,,,]

## Ditto location C
PCT.mean <- apply(P[,,,,3],c(2,3,4),mean)
PCT.min <- apply(P[,,,,3],c(2,3,4),rethinking::HPDI)[1,,,]
PCT.max <- apply(P[,,,,3],c(2,3,4),rethinking::HPDI)[2,,,]

## Calculate raw time-specific probabilities to be at A, B, C:
PATraw <- array(dim=c(n_species,n_treat,n_time))
PBTraw <- array(dim=c(n_species,n_treat,n_time))
PCTraw <- array(dim=c(n_species,n_treat,n_time))

## treat (infected, test) 0:0,0:1,1:0,1:1

d$treat <-   as.integer(ifelse(d$test==0 & d$infected==0,1,
                               ifelse(d$test==1 & d$infected==0,2,
                                      ifelse(d$test==0 & d$infected==1,3,4))))

sample_sizes <- table(d$species,d$treat)

for (i in 1:n_time) PATraw[,,i]<-tapply(as.data.frame(d)[,4+i],list(d$species,d$treat),function(x) sum(x=="A"))/sample_sizes
for (i in 1:n_time) PBTraw[,,i]<-tapply(as.data.frame(d)[,4+i],list(d$species,d$treat),function(x) sum(x=="B"))/sample_sizes
for (i in 1:n_time) PCTraw[,,i]<-tapply(as.data.frame(d)[,4+i],list(d$species,d$treat),function(x) sum(x=="C"))/sample_sizes


species <- c("Species EB","Species GF","Species GP")
treats <- c("UC","UT","IC","IT")

dg.time <- expand.grid(Species = species,
                       Treatment = treats,
                       Time = 1:10,
                       Area = c("A","B","C"))

dg.time <- transform(dg.time,
                     y = c(as.vector(PAT.mean),
                           as.vector(PBT.mean),
                           as.vector(PCT.mean)))

dg.time <- transform(dg.time,
                     ymin = c(as.vector(PAT.min),
                              as.vector(PBT.min),
                              as.vector(PCT.min)))

dg.time <- transform(dg.time,
                     ymax = c(as.vector(PAT.max),
                              as.vector(PBT.max),
                              as.vector(PCT.max)))

dg.time <- transform(dg.time,
                     yraw = c(as.vector(PATraw),
                              as.vector(PBTraw),
                              as.vector(PCTraw)))

ggplot(dg.time,aes(x=Time,y=y,color=Area)) +
  facet_grid(Species~Treatment) +
  geom_line() +
  geom_point(aes(x=Time,y=yraw),size=1) +
  geom_ribbon(aes(ymin=ymin,ymax=ymax,fill=Area),alpha=0.2,color=NA) +
  background_grid(major = "xy", minor = "xy") +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,5))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.5))+
  labs(x="\nTime (s)",y="Proportion present\n")


#### Figure 3: total proportions of time spent per area:

## Sorted by location, treatment, species
PX_mean <- as.vector(apply(Q,c(2,3,4),mean))
PX_min <- as.vector(apply(Q,c(2,3,4),rethinking::HPDI)[1,,,])
PX_max <- as.vector(apply(Q,c(2,3,4),rethinking::HPDI)[2,,,])


species <- c("Species EB","Species GF","Species GP")
treats <- c("UC","UT","IC","IT")

dg.loc <- expand.grid(Species = species,
                      Treatment = treats,
                      Area = c("A","B","C"))

dg.loc <- transform(dg.loc,
                    y = PX_mean,
                    ymin = PX_min,
                    ymax = PX_max)
dg.loc$TC_group <- as.factor(c(1,2,3,1,4,5,6,7,8,9,10,8,
                               11,12,13,14,15,13,16,17,18,16,17,18,
                               19,20,21,22,23,24,25,26,27,28,29,27))

sp <- levels(dg.loc$Species)
tr <- levels(dg.loc$Treatment)

ann_text <- data.frame(Species=c(sp[2],sp[2],sp[3],sp[3],sp[3],sp[3]),
                       Treatment=c(tr[2],tr[3],tr[1],tr[2],tr[3],tr[4]),
                       Area=c("A","B","B","A","A","A"),
                       y=rep(-0.05,6),
                       label=c("***","**","**","***","*","*"))

ggplot(dg.loc,aes(x=Treatment,y=y,color=Area)) +
  facet_grid(~Species) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=ymin,ymax=ymax),width=0.15) +
  background_grid(major = "xy", minor = "xy") +
  scale_y_continuous(limits=c(-0.07,1),breaks=seq(0,1,0.5))+
  labs(y="Proportion of time present\n") +
  geom_text(data = ann_text,label=ann_text$label,show.legend = FALSE) +
  geom_line(aes(group=TC_group))


### posterior probs of Pr(A)-Pr(B)>0

PAminB_mean <- apply(Q[,,,1]-Q[,,,2]>0,c(2,3),sum)/n_sim

#        [,1]   [,2]   [,3]   [,4]
# [1,] 0.8284 0.4326 0.1603 0.7800
# [2,] 0.0922 0.9994 0.0057 0.9244
# [3,] 0.0029 1.0000 0.9897 0.9624

## These "p-values" were used to anotate "significance" stars at the
## bottom of Figure 3
## *   = p<0.05 or p>0.95
## **  = p<0.01 or p>0.99
## *** = p<0.001 or p>0.999


### Posterior probs of the 12 T - C differences for A and B

PUTminUC <- apply(Q[,,2,]-Q[,,1,]>0,c(2,3),sum)/n_sim
PUTminUC # rows: species, cols: areas
#        [,1]   [,2]   [,3]
# [1,] 0.8804 0.9855 0.0094
# [2,] 1.0000 0.9569 0.0000
# [3,] 1.0000 0.0617 0.0002

PITminIC <- apply(Q[,,4,]-Q[,,3,]>0,c(2,3),sum)/n_sim
PITminIC # rows: species, cols: areas
#        [,1]   [,2]   [,3]
# [1,] 0.9805 0.7282 0.0270
# [2,] 0.9999 0.4087 0.0184
# [3,] 0.7855 0.7103 0.1573

## The NS differences were encoded as lines in Figure 3
## by giving NS T and C combinations the same level in
## a special factor TC_group. Lines are then drawn with
## geom_lines(aes(group=TC_group))

## Posterior probability that Pr(A)_UT - Pr(A)_UC > Pr(A)_IT - Pr(A)_IC

apply(Q[,,2,]-Q[,,1,]>Q[,,4,]-Q[,,3,],c(2,3),sum)/n_sim
