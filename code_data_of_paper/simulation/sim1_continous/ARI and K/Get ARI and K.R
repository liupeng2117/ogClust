library(mclust)
setwd("D:/research/Disease subtyping/simulation/compare with clustering/result/groups and genes/")
n=600
#----- skmeans expression 3 gamma 1 delta 3 -----#
load("skmeans_cd100_exp3_gamma1_delta3.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3
cl.skm3<-c()
ARI.skm3<-c()
ARI2.skm3<-c()
ARI3.skm3<-c()
for(i in 1:100){
  cl.skm3[i]<-res.skm[[i]][[1]] #group assign
  ARI.skm3[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign)
  ARI2.skm3[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign2)
  ARI3.skm3[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign3)
}
table(cl.skm3)
# cl.skm3
# 2  3 
# 1 99

#ARI
summary(ARI.skm3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2749  0.5706  0.5706  0.5676  0.5706  0.5706  
summary(ARI2.skm3)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0005828 -0.0005828 -0.0005828 -0.0005700 -0.0005828  0.0006949 
summary(ARI3.skm3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5706  1.0000  1.0000  0.9957  1.0000  1.0000 

# Genes: FNs FPs
TPs.skm3<-c()
TNs.skm3<-c()
FNs.skm3<-c()
FPs.skm3<-c()
for(i in 1:100){
  TPs.skm3[i]<-sum(which(res.skm[[i]][[3]]>=0.05) %in% 1:10)
  TNs.skm3[i]<-sum(which(res.skm[[i]][[3]]<0.05) %in% 11:100)
  FNs.skm3[i]<-sum(which(res.skm[[i]][[3]]<0.05) %in% 1:10)
  FPs.skm3[i]<-sum(which(res.skm[[i]][[3]]>=0.05) %in% 11:100)
}

summary(FNs.skm3)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   10      10      10      10      10      10 
summary(FPs.skm3)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  15      15      15      15      15      15 



#----- skmeans expression 3 gamma 0.5 delta 3 -----#
load("skmeans_cd100_exp3_gamma0.5_delta3.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3
cl.skm3<-c()
ARI.skm3<-c()
ARI2.skm3<-c()
ARI3.skm3<-c()
for(i in 1:100){
  cl.skm3[i]<-res.skm[[i]][[1]] #group assign
  ARI.skm3[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign)
  ARI2.skm3[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign2)
  ARI3.skm3[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign3)
}
table(cl.skm3)
# cl.skm3
# 2  3 
# 1 99

#ARI
summary(ARI.skm3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2749  0.5706  0.5706  0.5676  0.5706  0.5706  
summary(ARI2.skm3)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0005828 -0.0005828 -0.0005828 -0.0005700 -0.0005828  0.0006949 
summary(ARI3.skm3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5706  1.0000  1.0000  0.9957  1.0000  1.0000 

# Genes: FNs FPs
TPs.skm3<-c()
TNs.skm3<-c()
FNs.skm3<-c()
FPs.skm3<-c()
for(i in 1:100){
  TPs.skm3[i]<-sum(which(res.skm[[i]][[3]]>=0.05) %in% 1:10)
  TNs.skm3[i]<-sum(which(res.skm[[i]][[3]]<0.05) %in% 11:100)
  FNs.skm3[i]<-sum(which(res.skm[[i]][[3]]<0.05) %in% 1:10)
  FPs.skm3[i]<-sum(which(res.skm[[i]][[3]]>=0.05) %in% 11:100)
}

summary(FNs.skm3)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   10      10      10      10      10      10 
summary(FPs.skm3)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  15      15      15      15      15      15 

#----- skmeans expression 3 gamma 1 delta 1 -----#
load("skmeans_cd100_exp3_gamma1_delta1.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3
cl.skm3<-c()
ARI.skm3<-c()
ARI2.skm3<-c()
ARI3.skm3<-c()
for(i in 1:100){
  cl.skm3[i]<-res.skm[[i]][[1]] #group assign
  ARI.skm3[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign)
  ARI2.skm3[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign2)
  ARI3.skm3[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign3)
}
table(cl.skm3)
# cl.skm3
# 2  3 
# 1 99

#ARI
summary(ARI.skm3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2749  0.5706  0.5706  0.5676  0.5706  0.5706  
summary(ARI2.skm3)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0005828 -0.0005828 -0.0005828 -0.0005700 -0.0005828  0.0006949 
summary(ARI3.skm3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5706  1.0000  1.0000  0.9957  1.0000  1.0000 

# Genes: FNs FPs
TPs.skm3<-c()
TNs.skm3<-c()
FNs.skm3<-c()
FPs.skm3<-c()
for(i in 1:100){
  TPs.skm3[i]<-sum(which(res.skm[[i]][[3]]>=0.05) %in% 1:10)
  TNs.skm3[i]<-sum(which(res.skm[[i]][[3]]<0.05) %in% 11:100)
  FNs.skm3[i]<-sum(which(res.skm[[i]][[3]]<0.05) %in% 1:10)
  FPs.skm3[i]<-sum(which(res.skm[[i]][[3]]>=0.05) %in% 11:100)
}

summary(FNs.skm3)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   10      10      10      10      10      10 
summary(FPs.skm3)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  15      15      15      15      15      15 

#----- skmeans expression 3 gamma 0.5 delta 1 -----#
load("skmeans_cd100_exp3_gamma0.5_delta1.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3
cl.skm3<-c()
ARI.skm3<-c()
ARI2.skm3<-c()
ARI3.skm3<-c()
for(i in 1:100){
  cl.skm3[i]<-res.skm[[i]][[1]] #group assign
  ARI.skm3[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign)
  ARI2.skm3[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign2)
  ARI3.skm3[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign3)
}
table(cl.skm3)
# cl.skm3
# 2  3 
# 1 99

#ARI
summary(ARI.skm3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2749  0.5706  0.5706  0.5676  0.5706  0.5706  
summary(ARI2.skm3)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0005828 -0.0005828 -0.0005828 -0.0005700 -0.0005828  0.0006949 
summary(ARI3.skm3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5706  1.0000  1.0000  0.9957  1.0000  1.0000 

# Genes: FNs FPs
TPs.skm3<-c()
TNs.skm3<-c()
FNs.skm3<-c()
FPs.skm3<-c()
for(i in 1:100){
  TPs.skm3[i]<-sum(which(res.skm[[i]][[3]]>=0.05) %in% 1:10)
  TNs.skm3[i]<-sum(which(res.skm[[i]][[3]]<0.05) %in% 11:100)
  FNs.skm3[i]<-sum(which(res.skm[[i]][[3]]<0.05) %in% 1:10)
  FPs.skm3[i]<-sum(which(res.skm[[i]][[3]]>=0.05) %in% 11:100)
}

summary(FNs.skm3)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   10      10      10      10      10      10 
summary(FPs.skm3)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  15      15      15      15      15      15 

#----- skmeans expression 1 gamma 1 delta 3 -----#
load("skmeans_cd100_exp1_gamma1_delta3.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.skm1<-c()
ARI.skm1<-c()
ARI2.skm1<-c()
ARI3.skm1<-c()
for(i in 1:100){
  cl.skm1[i]<-res.skm[[i]][[1]] #group assign
  ARI.skm1[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign)
  ARI2.skm1[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign2)
  ARI3.skm1[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign3)
}
table(cl.skm1)
#  cl.skm3
# 2  3  6 
#60  6 34 

#ARI
summary(ARI.skm1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2185  0.2444  0.2603  0.3685  0.5617  0.6720 
summary(ARI2.skm1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2057  0.2543  0.7192  0.5601  0.7623  0.8217 
summary(ARI3.skm1)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0020946 -0.0007885  0.0014145  0.1342340  0.3501258  0.4175072 

# Genes: FNs, FPs
TPs.skm1<-c()
FNs.skm1<-c()
FNs.skm1<-c()
FPs.skm1<-c()
for(i in 1:100){
  TPs.skm1[i]<-sum(which(res.skm[[i]][[3]]>=0.05) %in% 1:10)
  TNs.skm1[i]<-sum(which(res.skm[[i]][[3]]<0.05) %in% 11:100)
  FNs.skm1[i]<-sum(which(res.skm[[i]][[3]]<0.05) %in% 1:10)
  FPs.skm1[i]<-sum(which(res.skm[[i]][[3]]>=0.05) %in% 11:100)
}

summary(FNs.skm1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0       0       0       0       0       0 
summary(FPs.skm1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00    0.00    0.00    5.43   15.00   15.00 


#----- skmeans expression 0.5 gamma 1 delta 3 -----#
load("skmeans_cd100_exp0.5_gamma1_delta3.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.skm0.5<-c()
ARI.skm0.5<-c()
ARI2.skm0.5<-c()
ARI3.skm0.5<-c()
for(i in 1:100){
  cl.skm0.5[i]<-res.skm[[i]][[1]] #group assign
  ARI.skm0.5[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign)
  ARI2.skm0.5[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign2)
  ARI3.skm0.5[i]<-adjustedRandIndex(res.skm[[i]][[4]],grp_assign3)
}
table(cl.skm0.5)
#  cl.skm3
# 2  3  6 
#75  0 25 

#ARI
summary(ARI.skm0.5)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1274  0.2058  0.2438  0.2225  0.2557  0.2752 
summary(ARI2.skm0.5)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2074  0.5641  0.7334  0.6235  0.7696  0.8278 
summary(ARI3.skm0.5)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0035533 -0.0015831 -0.0007359 -0.0004589  0.0002880  0.0059293 

# Genes: FNs, FPs
TPs.skm0.5<-c()
TNs.skm0.5<-c()
FNs.skm0.5<-c()
FPs.skm0.5<-c()
for(i in 1:100){
  TPs.skm0.5[i]<-sum(which(res.skm[[i]][[3]]>=0.05) %in% 1:10)
  TNs.skm0.5[i]<-sum(which(res.skm[[i]][[3]]<0.05) %in% 11:100)
  FNs.skm0.5[i]<-sum(which(res.skm[[i]][[3]]<0.05) %in% 1:10)
  FPs.skm0.5[i]<-sum(which(res.skm[[i]][[3]]>=0.05) %in% 11:100)
}

summary(FNs.skm0.5)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0       0       0       0       0       0 
summary(FPs.skm0.5)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0       0       0       0       0       0 

#----- pmbc expression 3 gamma 1 delta 3 -----#
load("pmbc_cd100_exp3_gamma1_delta3.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.pmbc3<-c()
ARI.pmbc3<-c()
ARI2.pmbc3<-c()
ARI3.pmbc3<-c()
for(i in 1:100){
  cl.pmbc3[i]<-res.pmbc[[i]][[1]] #group assign
  ARI.pmbc3[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign)
  ARI2.pmbc3[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign2)
  ARI3.pmbc3[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign3)
}
table(cl.pmbc3)
#  cl.pmbc3
#  3  6 
# 91  9 

#ARI
summary(ARI.pmbc3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2530  0.5706  0.5706  0.5495  0.5706  0.5706 
summary(ARI2.pmbc3)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0010856 -0.0005828 -0.0005828  0.0312585 -0.0005828  0.8339006 
summary(ARI3.pmbc3)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.002037  1.000000  1.000000  0.947125  1.000000  1.000000 

# Genes: FNs, FPs
TPs.pmbc3<-c()
TNs.pmbc3<-c()
FNs.pmbc3<-c()
FPs.pmbc3<-c()
for(i in 1:100){
  TPs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)>=0.05) %in% 1:10)
  TNs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)<0.05) %in% 11:100)
  FNs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)<0.05) %in% 1:10)
  FPs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)>=0.05) %in% 11:100)
}

summary(FNs.pmbc3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0    10.0    10.0     9.6    10.0    10.0 
summary(FPs.pmbc3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   15.00   15.00   14.43   15.00   16.00 

#----- pmbc expression 3 gamma 0.5 delta 3-----#
load("pmbc_cd100_exp3_gamma0.5_delta3.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.pmbc3<-c()
ARI.pmbc3<-c()
ARI2.pmbc3<-c()
ARI3.pmbc3<-c()
for(i in 1:100){
  cl.pmbc3[i]<-res.pmbc[[i]][[1]] #group assign
  ARI.pmbc3[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign)
  ARI2.pmbc3[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign2)
  ARI3.pmbc3[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign3)
}
table(cl.pmbc3)
#  cl.pmbc3
#  3  6 
# 91  9 

#ARI
summary(ARI.pmbc3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2530  0.5706  0.5706  0.5495  0.5706  0.5706 
summary(ARI2.pmbc3)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0010856 -0.0005828 -0.0005828  0.0312585 -0.0005828  0.8339006 
summary(ARI3.pmbc3)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.002037  1.000000  1.000000  0.947125  1.000000  1.000000 

# Genes: FNs, FPs
TPs.pmbc3<-c()
TNs.pmbc3<-c()
FNs.pmbc3<-c()
FPs.pmbc3<-c()
for(i in 1:100){
  TPs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)>=0.05) %in% 1:10)
  TNs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)<0.05) %in% 11:100)
  FNs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)<0.05) %in% 1:10)
  FPs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)>=0.05) %in% 11:100)
}

summary(FNs.pmbc3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0    10.0    10.0     9.6    10.0    10.0 
summary(FPs.pmbc3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   15.00   15.00   14.43   15.00   16.00 

#----- pmbc expression 3 gamma 1 delta 1-----#
load("pmbc_cd100_exp3_gamma1_delta1.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.pmbc3<-c()
ARI.pmbc3<-c()
ARI2.pmbc3<-c()
ARI3.pmbc3<-c()
for(i in 1:100){
  cl.pmbc3[i]<-res.pmbc[[i]][[1]] #group assign
  ARI.pmbc3[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign)
  ARI2.pmbc3[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign2)
  ARI3.pmbc3[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign3)
}
table(cl.pmbc3)
#  cl.pmbc3
#  3  6 
# 91  9 

#ARI
summary(ARI.pmbc3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2530  0.5706  0.5706  0.5495  0.5706  0.5706 
summary(ARI2.pmbc3)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0010856 -0.0005828 -0.0005828  0.0312585 -0.0005828  0.8339006 
summary(ARI3.pmbc3)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.002037  1.000000  1.000000  0.947125  1.000000  1.000000 

# Genes: FNs, FPs
TPs.pmbc3<-c()
TNs.pmbc3<-c()
FNs.pmbc3<-c()
FPs.pmbc3<-c()
for(i in 1:100){
  TPs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)>=0.05) %in% 1:10)
  TNs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)<0.05) %in% 11:100)
  FNs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)<0.05) %in% 1:10)
  FPs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)>=0.05) %in% 11:100)
}

summary(FNs.pmbc3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0    10.0    10.0     9.6    10.0    10.0 
summary(FPs.pmbc3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   15.00   15.00   14.43   15.00   16.00 

#----- pmbc expression 3 gamma 0.5 delta 1-----#
load("pmbc_cd100_exp3_gamma0.5_delta1.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.pmbc3<-c()
ARI.pmbc3<-c()
ARI2.pmbc3<-c()
ARI3.pmbc3<-c()
for(i in 1:100){
  cl.pmbc3[i]<-res.pmbc[[i]][[1]] #group assign
  ARI.pmbc3[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign)
  ARI2.pmbc3[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign2)
  ARI3.pmbc3[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign3)
}
table(cl.pmbc3)
#  cl.pmbc3
#  3  6 
# 91  9 

#ARI
summary(ARI.pmbc3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2530  0.5706  0.5706  0.5495  0.5706  0.5706 
summary(ARI2.pmbc3)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0010856 -0.0005828 -0.0005828  0.0312585 -0.0005828  0.8339006 
summary(ARI3.pmbc3)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.002037  1.000000  1.000000  0.947125  1.000000  1.000000 

# Genes: FNs, FPs
TPs.pmbc3<-c()
TNs.pmbc3<-c()
FNs.pmbc3<-c()
FPs.pmbc3<-c()
for(i in 1:100){
  TPs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)>=0.05) %in% 1:10)
  TNs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)<0.05) %in% 11:100)
  FNs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)<0.05) %in% 1:10)
  FPs.pmbc3[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)>=0.05) %in% 11:100)
}

summary(FNs.pmbc3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0    10.0    10.0     9.6    10.0    10.0 
summary(FPs.pmbc3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   15.00   15.00   14.43   15.00   16.00 

#----- pmbc expression 1 -----#
load("pmbc_cd100_exp1_gamma1_delta3.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.pmbc1<-c()
ARI.pmbc1<-c()
ARI2.pmbc1<-c()
ARI3.pmbc1<-c()
for(i in 1:100){
  cl.pmbc1[i]<-res.pmbc[[i]][[1]] #group assign
  ARI.pmbc1[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign)
  ARI2.pmbc1[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign2)
  ARI3.pmbc1[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign3)
}
table(cl.pmbc1)
#  cl.pmbc1
# 2  3 
# 87 13 

#ARI
summary(ARI.pmbc1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1541  0.1804  0.2500  0.2317  0.2625  0.4055 
summary(ARI2.pmbc1)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0020850 -0.0002835  0.7478069  0.4967352  0.7917523  0.8461434
summary(ARI3.pmbc1)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0021249 -0.0011006  0.0004893  0.1331263  0.3371496  0.7134314 

# Genes: FNs, FPs
TPs.pmbc1<-c()
TNs.pmbc1<-c()
FNs.pmbc1<-c()
FPs.pmbc1<-c()
for(i in 1:100){
  TPs.pmbc1[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)>=0.05) %in% 1:10)
  TNs.pmbc1[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)<0.05) %in% 11:100)
  FNs.pmbc1[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)<0.05) %in% 1:10)
  FPs.pmbc1[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)>=0.05) %in% 11:100)
}

summary(FNs.pmbc1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0     0.0     0.0     3.6    10.0    10.0 
summary(FPs.pmbc1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    0.00    4.25   10.00   15.00 

#----- pmbc expression 0.5 -----#
load("pmbc_cd100_exp0.5_gamma1_delta3.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.pmbc0.5<-c()
ARI.pmbc0.5<-c()
ARI2.pmbc0.5<-c()
ARI3.pmbc0.5<-c()
for(i in 1:100){
  cl.pmbc0.5[i]<-res.pmbc[[i]][[1]] #group assign
  ARI.pmbc0.5[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign)
  ARI2.pmbc0.5[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign2)
  ARI3.pmbc0.5[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),grp_assign3)
}
table(cl.pmbc0.5)
#  cl.pmbc0.5
# 2
# 100

#ARI
summary(ARI.pmbc0.5)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.2385  0.2539  0.1990  0.2636  0.2818
summary(ARI2.pmbc0.5)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.7206  0.7652  0.5992  0.7932  0.8461
summary(ARI3.pmbc0.5)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0021249 -0.0012907 -0.0003649 -0.0004469  0.0000000  0.0032930   

# Genes: FNs, FPs
TPs.pmbc0.5<-c()
TNs.pmbc0.5<-c()
FNs.pmbc0.5<-c()
FPs.pmbc0.5<-c()
for(i in 1:100){
  TPs.pmbc0.5[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)>=0.05) %in% 1:10)
  TNs.pmbc0.5[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)<0.05) %in% 11:100)
  FNs.pmbc0.5[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)<0.05) %in% 1:10)
  FPs.pmbc0.5[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)>=0.05) %in% 11:100)
}

summary(FNs.pmbc0.5)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0     0.0     0.0     2.3     0.0    10.0 
summary(FPs.pmbc0.5)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.00    0.00    0.00    0.18    0.00    2.00 

#----- gfmr expression 3 gamma 1 delta 3 -----#
load("gfmr_cd100_exp3_gamma1_delta3.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.gfmr3<-c()
ARI.gfmr3<-c()
ARI2.gfmr3<-c()
ARI3.gfmr3<-c()
for(i in 1:100){
  cl.gfmr3[i]<-res.gfmr[[i]][[1]] #group assign
  ARI.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign)
  ARI2.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign2)
  ARI3.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign3)
}
table(cl.gfmr3)
# cl.gfmr3
# 2  3 
# 99  1 

#ARI
summary(ARI.gfmr3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1814  0.2098  0.2176  0.2182  0.2264  0.2512 
summary(ARI2.gfmr3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5468  0.6341  0.6555  0.6569  0.6828  0.7565 
summary(ARI3.gfmr3)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0021852 -0.0014587 -0.0008498 -0.0005356  0.0001007  0.0056236 

# Genes: FNs, FPs
TPs.gfmr3<-c()
TNs.gfmr3<-c()
FNs.gfmr3<-c()
FPs.gfmr3<-c()
for(i in 1:100){
  TPs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])>=0.05) %in% 1:10)
  TNs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])<0.05) %in% 11:100)
  FNs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])<0.05) %in% 1:10)
  FPs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])>=0.05) %in% 11:100)
}

summary(FNs.gfmr3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    0.00    0.03    0.00    3.00 
summary(FPs.gfmr3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    0.00    0.92    1.00    7.00 

#----- gfmr expression 3 gamma 1 delta 3 -----#
load("gfmr_cd100_exp3_gamma1_delta3_new.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.gfmr3<-c()
ARI.gfmr3<-c()
ARI2.gfmr3<-c()
ARI3.gfmr3<-c()
for(i in 1:100){
  cl.gfmr3[i]<-res.gfmr[[i]][[1]] #group assign
  ARI.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign)
  ARI2.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign2)
  ARI3.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign3)
}
table(cl.gfmr3)
# cl.gfmr3
# 2  3 
# 99  1 

#ARI
summary(ARI.gfmr3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2191  0.2467  0.2530  0.2524  0.2615  0.2866 
summary(ARI2.gfmr3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6609  0.7449  0.7623  0.7613  0.7888  0.8647 
summary(ARI3.gfmr3)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.0021699 -0.0016799 -0.0009908 -0.0009475 -0.0003805  0.0013664

# Genes: FNs, FPs
TPs.gfmr3<-c()
TNs.gfmr3<-c()
FNs.gfmr3<-c()
FPs.gfmr3<-c()
for(i in 1:100){
  TPs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])>=0.05) %in% 1:10)
  TNs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])<0.05) %in% 11:100)
  FNs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])<0.05) %in% 1:10)
  FPs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])>=0.05) %in% 11:100)
}

summary(FNs.gfmr3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       0       0       0       0       0
summary(FPs.gfmr3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    1.00    1.30    2.25    4.00
#----- gfmr expression 3 gamma 0.5 delta 3 -----#
load("gfmr_cd100_exp3_gamma0.5_delta3.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.gfmr3<-c()
ARI.gfmr3<-c()
ARI2.gfmr3<-c()
ARI3.gfmr3<-c()
for(i in 1:100){
  cl.gfmr3[i]<-res.gfmr[[i]][[1]] #group assign
  ARI.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign)
  ARI2.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign2)
  ARI3.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign3)
}
table(cl.gfmr3)
#cl.gfmr3
#2 
#100 

#ARI
summary(ARI.gfmr3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1814  0.2098  0.2176  0.2182  0.2264  0.2512 
summary(ARI2.gfmr3)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5468  0.6341  0.6555  0.6569  0.6828  0.7565 
summary(ARI3.gfmr3)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0021852 -0.0014587 -0.0008498 -0.0005356  0.0001007  0.0056236 

# Genes: FNs, FPs
TPs.gfmr3<-c()
TNs.gfmr3<-c()
FNs.gfmr3<-c()
FPs.gfmr3<-c()
for(i in 1:100){
  TPs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])>=0.05) %in% 1:10)
  TNs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])<0.05) %in% 11:100)
  FNs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])<0.05) %in% 1:10)
  FPs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])>=0.05) %in% 11:100)
}

summary(FNs.gfmr3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    0.00    0.03    0.00    3.00 
summary(FPs.gfmr3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    0.00    0.92    1.00    7.00 

#----- gfmr expression 3 gamma 1 delta 1 -----#
load("gfmr_cd100_exp3_gamma1_delta1.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.gfmr3<-c()
ARI.gfmr3<-c()
ARI2.gfmr3<-c()
ARI3.gfmr3<-c()
for(i in 1:100){
  cl.gfmr3[i]<-res.gfmr[[i]][[1]] #group assign
  ARI.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign)
  ARI2.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign2)
  ARI3.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign3)
}
table(cl.gfmr3)
#cl.gfmr3
#2 
#100 

#ARI
summary(ARI.gfmr3[FNs.gfmr3<9])
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.003926 0.147474 0.171808 0.160489 0.202373 0.264771
summary(ARI2.gfmr3[FNs.gfmr3<9])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01822 0.44803 0.52727 0.49341 0.61297 0.79771
summary(ARI3.gfmr3[FNs.gfmr3<9])
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0021986 -0.0014857 -0.0010259 -0.0003985  0.0002327  0.0102941

# Genes: FNs, FPs
TPs.gfmr3<-c()
TNs.gfmr3<-c()
FNs.gfmr3<-c()
FPs.gfmr3<-c()
for(i in 1:100){
  TPs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])>=0.05) %in% 1:10)
  TNs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])<0.05) %in% 11:100)
  FNs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])<0.05) %in% 1:10)
  FPs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])>=0.05) %in% 11:100)
}

summary(FNs.gfmr3[FNs.gfmr3<9])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   3.000   3.000   3.648   4.000   8.000  
summary(FPs.gfmr3[FNs.gfmr3<9])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.0000  0.0000  0.1538  0.0000  2.0000 

#----- gfmr expression 3 gamma 1 delta 1 -----#
load("gfmr_cd100_exp3_gamma0.5_delta1.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.gfmr3<-c()
ARI.gfmr3<-c()
ARI2.gfmr3<-c()
ARI3.gfmr3<-c()
for(i in 1:100){
  cl.gfmr3[i]<-res.gfmr[[i]][[1]] #group assign
  ARI.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign)
  ARI2.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign2)
  ARI3.gfmr3[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign3)
}
table(cl.gfmr3)
#cl.gfmr3
#2 
#100 

#ARI
summary(ARI.gfmr3[FNs.gfmr3<9])
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000324 0.087221 0.153845 0.135781 0.192792 0.250514
summary(ARI2.gfmr3[FNs.gfmr3<9])
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.001715 0.288754 0.463781 0.419006 0.578176 0.756495 
summary(ARI3.gfmr3[FNs.gfmr3<9])
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-2.095e-03 -1.650e-03 -1.021e-03  1.483e-04  3.984e-05  2.466e-02  

# Genes: FNs, FPs
TPs.gfmr3<-c()
TNs.gfmr3<-c()
FNs.gfmr3<-c()
FPs.gfmr3<-c()
for(i in 1:100){
  TPs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])>=0.05) %in% 1:10)
  TNs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])<0.05) %in% 11:100)
  FNs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])<0.05) %in% 1:10)
  FPs.gfmr3[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])>=0.05) %in% 11:100)
}

summary(FNs.gfmr3[FNs.gfmr3<9])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    0.00    0.03    0.00    3.00 
summary(FPs.gfmr3[FNs.gfmr3<9])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.0000  0.0000  0.2083  0.0000  2.0000 

#----- gfmr expression 1 -----#
load("gfmr_cd100_exp1_gamma1_delta3.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.gfmr1<-c()
ARI.gfmr1<-c()
ARI2.gfmr1<-c()
ARI3.gfmr1<-c()
for(i in 1:100){
  cl.gfmr1[i]<-res.gfmr[[i]][[1]] #group assign
  ARI.gfmr1[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign)
  ARI2.gfmr1[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign2)
  ARI3.gfmr1[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign3)
}
table(cl.gfmr1)
# cl.gfmr1
# 2  3 
# 99  1 

#ARI
summary(ARI.gfmr1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1814  0.2098  0.2176  0.2179  0.2263  0.2512 
summary(ARI2.gfmr1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5468  0.6341  0.6555  0.6569  0.6828  0.7565 
summary(ARI3.gfmr1)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0021269 -0.0014652 -0.0008634 -0.0005461  0.0001007  0.0056236 

# Genes: FNs, FPs
TPs.gfmr1<-c()
TNs.gfmr1<-c()
FNs.gfmr1<-c()
FPs.gfmr1<-c()
for(i in 1:100){
  TPs.gfmr1[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])>=0.05) %in% 1:10)
  TNs.gfmr1[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])<0.05) %in% 11:100)
  FNs.gfmr1[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])<0.05) %in% 1:10)
  FPs.gfmr1[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])>=0.05) %in% 11:100)
}

summary(FNs.gfmr1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    0.00    0.04    0.00    4.00 
summary(FPs.gfmr1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    0.00    0.68    1.00    6.00 
#----- gfmr expression 0.5 -----#
load("gfmr_cd100_exp0.5_gamma1_delta3.rdata")
grp_assign2=rep(1,n)
grp_assign2[grp_assign %in% c(2,4,6)]<-2
grp_assign3=rep(1,n)
grp_assign3[grp_assign %in% c(3,4)]<-2
grp_assign3[grp_assign %in% c(5,6)]<-3

cl.gfmr0.5<-c()
ARI.gfmr0.5<-c()
ARI2.gfmr0.5<-c()
ARI3.gfmr0.5<-c()
for(i in 1:100){
  cl.gfmr0.5[i]<-res.gfmr[[i]][[1]] #group assign
  ARI.gfmr0.5[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign)
  ARI2.gfmr0.5[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign2)
  ARI3.gfmr0.5[i]<-adjustedRandIndex(apply(res.gfmr[[i]][[4]],1,which.max),grp_assign3)
}
table(cl.gfmr0.5)
# cl.gfmr0.5
# 2  3 
# 99  1 

#ARI
summary(ARI.gfmr0.5)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1798  0.2095  0.2176  0.2177  0.2263  0.2512 
summary(ARI2.gfmr0.5)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5388  0.6288  0.6555  0.6549  0.6828  0.7565 
summary(ARI3.gfmr0.5)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0021852 -0.0015010 -0.0008498 -0.0005584  0.0001007  0.0056236 

# Genes: FNs, FPs
TPs.gfmr0.5<-c()
TNs.gfmr0.5<-c()
FNs.gfmr0.5<-c()
FPs.gfmr0.5<-c()
for(i in 1:100){
  TPs.gfmr0.5[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])>=0.05) %in% 1:10)
  TNs.gfmr0.5[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])<0.05) %in% 11:100)
  FNs.gfmr0.5[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])<0.05) %in% 1:10)
  FPs.gfmr0.5[i]<-sum(which(abs(res.gfmr[[i]][[3]][4:103])>=0.05) %in% 11:100)
}

summary(FNs.gfmr0.5)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    0.00    0.04    0.00    4.00 
summary(FPs.gfmr0.5)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    0.00    0.74    1.00    6.00 


#----- scluster expression 3 -----#
load("scluster_cd100_exp3_gamma1_delta3.rdata")
cl.scluster3<-c()
ARI.scluster3<-c()
ARI2.scluster3<-c()
#ARI3.scluster3<-c()
for(i in 1:100){
  cl.scluster3[i]<-res.scluster[[i]][[1]] #group assign
  ARI.scluster3[i]<-adjustedRandIndex(res.scluster[[i]][[3]],grp_assign)
  ARI2.scluster3[i]<-adjustedRandIndex(res.scluster[[i]][[3]],grp_assign2)
  #ARI3.scluster3[i]<-adjustedRandIndex(res.scluster[[i]][[3]],grp_assign3)
}
table(cl.scluster3)
#cl.scluster3
#2  3 
#99  1 

#ARI
summary(ARI.scluster3)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.003001 -0.001528 -0.001167 -0.001000 -0.000592  0.000946 
summary(ARI2.scluster3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6664  0.7623  0.7799  0.7781  0.8037  0.8400  

# Genes: FNs FPs
TPs.scluster3<-c()
TNs.scluster3<-c()
FNs.scluster3<-c()
FPs.scluster3<-c()
for(i in 1:100){
  TPs.scluster3[i]<-sum(res.scluster[[i]][[2]] %in% 1:10)
  TNs.scluster3[i]<-sum((1:100)[-res.scluster[[i]][[2]]] %in% 11:100)
  FNs.scluster3[i]<-sum((1:100)[-res.scluster[[i]][[2]]] %in% 1:10)
  FPs.scluster3[i]<-sum(res.scluster[[i]][[2]] %in% 11:100)
}

summary(FNs.scluster3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    0.00    0.01    0.00    1.00 
summary(FPs.scluster3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    0.00    5.00    5.51   10.00   30.00 

#----- scluster expression 3 -----#
load("scluster_cd100_exp1_gamma1_delta3.rdata")
cl.scluster1<-c()
ARI.scluster1<-c()
ARI2.scluster1<-c()
#ARI3.scluster1<-c()
for(i in 1:100){
  cl.scluster1[i]<-res.scluster[[i]][[1]] #group assign
  ARI.scluster1[i]<-adjustedRandIndex(res.scluster[[i]][[3]],grp_assign)
  ARI2.scluster1[i]<-adjustedRandIndex(res.scluster[[i]][[3]],grp_assign2)
  #ARI3.scluster1[i]<-adjustedRandIndex(res.scluster[[i]][[3]],grp_assign3)
}
table(cl.scluster1)
#cl.scluster1
#2  3 
#99  1 

#ARI
summary(ARI.scluster1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1778  0.2012  0.2134  0.2153  0.2227  0.3470 
summary(ARI2.scluster1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6664  0.7565  0.7770  0.7769  0.7992  0.8461 

# Genes: FNs FPs
TPs.scluster1<-c()
TNs.scluster1<-c()
FNs.scluster1<-c()
FPs.scluster1<-c()
for(i in 1:100){
  TPs.scluster1[i]<-sum(res.scluster[[i]][[2]] %in% 1:10)
  TNs.scluster1[i]<-sum((1:100)[-res.scluster[[i]][[2]]] %in% 11:100)
  FNs.scluster1[i]<-sum((1:100)[-res.scluster[[i]][[2]]] %in% 1:10)
  FPs.scluster1[i]<-sum(res.scluster[[i]][[2]] %in% 11:100)
}

summary(FNs.scluster1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       0       0       0       0       0  
summary(FPs.scluster1)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0     0.0     5.0     8.1    10.0    30.0 

#----- scluster expression 3 -----#
load("scluster_cd100_exp0.5_gamma1_delta3.rdata")
cl.scluster0.5<-c()
ARI.scluster0.5<-c()
ARI2.scluster0.5<-c()
#ARI3.scluster0.5<-c()
for(i in 1:100){
  cl.scluster0.5[i]<-res.scluster[[i]][[1]] #group assign
  ARI.scluster0.5[i]<-adjustedRandIndex(res.scluster[[i]][[3]],grp_assign)
  ARI2.scluster0.5[i]<-adjustedRandIndex(res.scluster[[i]][[3]],grp_assign2)
  #ARI3.scluster0.5[i]<-adjustedRandIndex(res.scluster[[i]][[3]],grp_assign3)
}
table(cl.scluster0.5)
# cl.scluster0.5
# 2 
# 100 

#ARI
summary(ARI.scluster0.5)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1835  0.2042  0.2165  0.2207  0.2323  0.3470 
summary(ARI2.scluster0.5)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6664  0.7507  0.7799  0.7745  0.7977  0.8461  

# Genes: FNs FPs
TPs.scluster0.5<-c()
TNs.scluster0.5<-c()
FNs.scluster0.5<-c()
FPs.scluster0.5<-c()
for(i in 1:100){
  TPs.scluster0.5[i]<-sum(res.scluster[[i]][[2]] %in% 1:10)
  TNs.scluster0.5[i]<-sum((1:100)[-res.scluster[[i]][[2]]] %in% 11:100)
  FNs.scluster0.5[i]<-sum((1:100)[-res.scluster[[i]][[2]]] %in% 1:10)
  FPs.scluster0.5[i]<-sum(res.scluster[[i]][[2]] %in% 11:100)
}

summary(FNs.scluster0.5)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       0       0       0       0       0  
summary(FPs.scluster0.5)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       0      10      10      20      30  
