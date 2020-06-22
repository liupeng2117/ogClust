rm(list=ls())
source("~/subtyping/simulation/code/gen.data.K3.R")
# Simulation ------------------------------------------------
## Data generation --------------------
#============================
# 600 subjects, 2 subgroups, 6 clusters, 2 clinical variables x1, x2
# 100 genes G1 to G100
# G1 and G10 form two clusters and related to outcome
# G11 and G25 forms another three clusters but not related
# G26 to G100 are noise
# 10% outliers on Y
#============================
set.seed(123)
n.rep=100
NG=1000
n=600 # sample
np=2 # 2 covariates
beta1=1 # coefficient of clinical var1
beta2=1 # coefficient of clinical var2
gamma1=c(1,1,1,1,1,0,0,0,0,0,-1,-1,-1,-1,-1) # 1st set coefficient of genes 1 to 15, others are 0
gamma2=c(0,0,0,0,0,1,1,1,1,1,-1,-1,-1,-1,-1) # 2nd set coefficient of genes 1 to 15, others are 0
miu=c(1,4,7) # overall mean
sigma2=1 # overall variance
theta=list(beta1=beta1,
           beta2=beta2,
           gamma1=gamma1,
           gamma2=gamma2,
           miu=miu,
           sigma2=sigma2)

c1_index<-sample(n,n/3)
c2_index<-sample((1:n)[-c1_index],n/3)
c3_index<-(1:n)[-c(c1_index,c2_index)]

c4_index<-sample(n,n/3) #index for cluster 3
c5_index<-sample((1:n)[-c4_index],n/3) #index for cluster 4
c6_index<-(1:n)[c(-c4_index,-c5_index)] #index for cluster 5
cluster_index=list(
  c1_index=c1_index,
  c2_index=c2_index,
  c3_index=c3_index,
  c4_index=c4_index,
  c5_index=c5_index,
  c6_index=c6_index
)

gen.res<-gen.data.K3(n.rep, NG, n, theta, cluster_index)
data<-gen.res$data

library(doParallel) # foreach
cl<-makeCluster(10)
registerDoParallel(cl)

library(KLClust)
library(NbClust)
library(mclust)
res.skm <- foreach(rep = 1:n.rep, .packages=c("mclust","KLClust","NbClust")) %dopar% {
  print(paste("data",rep))
  df<-data[[rep]]
  G=as.matrix(df[,3:(NG+2)])
  Gs=sapply(1:ncol(G), function(i) (G[,i]-mean(G[,i]))/sd(G[,i]))
  #Y=df$Y
  #choose k
  optimal_k<-NbClust(data=Gs,min.nc = 2, max.nc = 9,method="kmeans",index="ch")$Best.nc["Number_clusters"]
  #kl<-KL.S4(t(Gs),lambda_list =  list(seq(3,5,0.5),seq(3,5,0.5),seq(3,5,0.5),seq(3,5,0.5)), 
  #          k_vector=c(2,3,6,9),num.cores = 1)
  #choose lambda
  optimal_lambda<-KMeansSparseCluster.permute(Gs,K=optimal_k,wbounds=seq(3,5,0.5))$bestw
  cltrain<-KMeansSparseCluster(Gs,K=optimal_k,wbounds=optimal_lambda)
  #cltrain<-KMeansSparseCluster(Gs,K=kl$optimal_k,wbounds=kl$optimal_lambda)
  weight<-cltrain[[1]]$ws
  cluster<-cltrain[[1]]$Cs
  return(list(optimal_k, optimal_lambda, weight, cluster))
}

save(res.skm,gen.res,file="skmeans_K3_1000genes_9groups_exp3_gamma1_delta3.rdata")

est.K<-c()
ARI<-c()
for(i in 1:10){
  est.K[i]<-res.skm[[i]][[1]]
  ARI[i]<-adjustedRandIndex(res.skm[[i]][[4]],gen.res$grp_assign[i,])
}
table(est.K)
ARI
median(ARI)

# Genes: FNs FPs
TPs<-c()
TNs<-c()
FNs<-c()
FPs<-c()
for(i in 1:10){
  TPs[i]<-sum(which(abs(res.skm[[i]][[3]])>=0.05) %in% 1:15)
  TNs[i]<-sum(which(abs(res.skm[[i]][[3]])<0.05) %in% 16:1000)
  FNs[i]<-sum(which(abs(res.skm[[i]][[3]])<0.05) %in% 1:15)
  FPs[i]<-sum(which(abs(res.skm[[i]][[3]])>=0.05) %in% 16:1000)
}

summary(FNs)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#15      15      15      15      15      15 
summary(FPs)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#15      15      15      15      15      15 
