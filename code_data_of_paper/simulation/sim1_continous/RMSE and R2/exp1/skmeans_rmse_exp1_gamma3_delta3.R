rm(list=ls())
source("~/subtyping/simulation/code/gen.data.K3.R")
#source("D:/research/Disease subtyping/simulation/gen.data.K3.R")
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
n.rep=10
NG=1000
n=600 # sample
np=2 # 2 covariates
beta1=1 # coefficient of clinical var1
beta2=1 # coefficient of clinical var2
gamma1=c(3,3,3,3,3,0,0,0,0,0,-3,-3,-3,-3,-3) # 1st set coefficient of genes 1 to 15, others are 0
gamma2=c(0,0,0,0,0,3,3,3,3,3,-3,-3,-3,-3,-3) # 2nd set coefficient of genes 1 to 15, others are 0
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

gen.res<-gen.data.K3(n.rep, NG, n, theta, cluster_index, exp=1)
data<-gen.res$data

library(doParallel) # foreach
cl<-makeCluster(10)
registerDoParallel(cl)

# Sparse kmeans
# cluster -----------------------------------------------------------------
library(sparcl)
library(caret)
#library(KLClust)
flds <- createFolds(1:n, k = 5, list = TRUE, returnTrain = FALSE)

#10-fold cross validation
res.skmeans<-foreach(rep = 1:n.rep, .packages=c("mclust","KLClust","NbClust")) %dopar% {
  print(rep)
  Y_prd=rep(0,n)
  grp_assign_prd=rep(0,n)
  for(f in 1:5){
    #train
    df.train<-data[[rep]][-flds[[f]],]
    #   T_train=df.train$T
    X1_train=df.train$X1
    X2_train=df.train$X2
    G_train=as.matrix(df.train[,3:(NG+2)])
    G_train=sapply(1:ncol(G_train), function(x) (G_train[,x]-mean(G_train[,x]))/sd(G_train[,x]))
    Y_train=df.train$Y
    
    #test
    df.test<-data[[rep]][flds[[f]],]
    #   T_test=df.test$T
    X1_test=df.test$X1
    X2_test=df.test$X2
    G_test=as.matrix(df.test[,3:(NG+2)])
    G_test=sapply(1:ncol(G_test), function(x) (G_test[,x]-mean(G_test[,x]))/sd(G_test[,x]))
    Y_test=df.test$Y
    
    #cltrain.per<-KMeansSparseCluster.permute(G_train,K=6)
    #kl<-KL.S4(t(G_train),lambda_list = list(seq(3,5,0.5),seq(3,5,0.5),seq(3,5,0.5),seq(3,5,0.5)), 
    #      k_vector=c(2,3,6,9),num.cores = 1)
    #cltrain<-KMeansSparseCluster(G_train,K=kl$optimal_k,wbounds=kl$optimal_lambda)
    optimal_k<-NbClust(data=G_train,min.nc = 2, max.nc = 9,method="kmeans",index="ch")$Best.nc["Number_clusters"]
    optimal_lambda<-KMeansSparseCluster.permute(G_train,K=optimal_k,wbounds=seq(3,5,0.5))$bestw
    cltrain<-KMeansSparseCluster(G_train,K=optimal_k,wbounds=optimal_lambda)
    
    weight<-cltrain[[1]]$ws
    grp_assign_train<-cltrain[[1]]$Cs
    n.k=length(unique(grp_assign_train))
    #distance matrix
    dist.matrix=matrix(,nrow=nrow(G_test),ncol=n.k)
    for(k in 1: n.k){
      miu_k<-apply(G_train[grp_assign_train==k,],2,mean)
      dist.matrix[,k]<-apply(G_test,1,function(x) weight %*% (x-miu_k)^2)
    }
    #group assignment of test data
    grp_assign_test=apply(dist.matrix, 1, which.min)
    grp_assign_prd[flds[[f]]]=grp_assign_test
    #calculate the y_est
    y_prd_cv=Y_test
    for(k in 1: n.k){
      y_k_train=Y_train[grp_assign_train==k]
      x1_k_train=X1_train[grp_assign_train==k]
      x2_k_train=X2_train[grp_assign_train==k]
      fit_k<-lm(y_k_train ~ x1_k_train + x2_k_train)
      y_k_test=Y_test[grp_assign_test==k]
      x1_k_test=X1_test[grp_assign_test==k]
      x2_k_test=X2_test[grp_assign_test==k]
      y_prd_cv[grp_assign_test==k]=fit_k$coefficients %*% rbind(1, x1_k_test, x2_k_test)
    }
    Y_prd[flds[[f]]]=y_prd_cv
  }
  #calculate RMSE
  RMSE=sqrt(sum((data[[rep]]$Y-Y_prd)^2)/n)
  R2=1-sum((data[[rep]]$Y-Y_prd)^2)/sum((data[[rep]]$Y-mean(data[[rep]]$Y))^2)
  return(list(RMSE=RMSE, R2=R2, grp_assign_prd=grp_assign_prd))
}

save(res.skmeans, gen.res, file="skmeans_rmse_K3_1000genes_9groups_exp1_gamma3_delta3.rdata")

library(mclust)
RMSE.all<-c()
R2<-c()
ARI.all<-c()
for(rep in 1:n.rep){
  RMSE.all[rep]<-res.skmeans[[rep]]$RMSE
  R2[rep]<-res.skmeans[[rep]]$R2
  grp_assign_prd<-res.skmeans[[rep]]$grp_assign_prd
  ARI.all[rep]<-adjustedRandIndex(gen.res$grp_assign[rep,], grp_assign_prd)
  print(table(gen.res$grp_assign[rep,], grp_assign_prd))
}
summary(RMSE.all)
summary(R2)
summary(ARI.all)