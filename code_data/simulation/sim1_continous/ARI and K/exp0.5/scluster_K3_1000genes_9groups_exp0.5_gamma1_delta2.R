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
miu=c(1,3,5) # overall mean
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

gen.res<-gen.data.K3(n.rep, NG, n, theta, cluster_index, exp=0.5)
data<-gen.res$data

library(doParallel) # foreach
cl<-makeCluster(10)
registerDoParallel(cl)

library(NbClust)
res.scluster <- foreach(rep = 1:n.rep, .packages=c("NbClust")) %dopar% {
  print(paste("data",rep))
  df<-data[[rep]]
  G=as.matrix(df[,3:(NG+2)])
  X1=df[,1]
  X2=df[,2]
  Y=df[,(NG+4)]
  Gs=sapply(1:ncol(G), function(i) (G[,i]-mean(G[,i]))/sd(G[,i]))
  
  # select the outcome related genes
  ps<-c()
  for(i in 1:ncol(Gs)){
    ps[i]<-summary(lm(Y ~ Gs[,i]))$coefficients[2,4]
  }
  ngs<-c(10,15,30,60,100)
  RMSE=c()
  # cross validation to select the best threshold
  
  for(o in 1:length(ngs)){
    ngsm<-ngs[o]
    Gs2<-Gs[,order(ps)[1:ngsm]]
    
    # 5 fold cross validation
    #parallel computing
    #library(doParallel) # foreach
    #cl<-makeCluster(4)
    #registerDoParallel(cl)
    
    # 5 fold
    n=nrow(Gs2)
    library(caret)
    set.seed(123)
    flds <- createFolds(1:n, k = 5, list = TRUE, returnTrain = FALSE)
    
    cvres.scluster<-list()
    for(f in 1:length(flds)){
      #train and test data
      Gs2_train<-Gs2[-flds[[f]],] #train g
      Gs2_test<-Gs2[flds[[f]],] #test g
      X1_train<-X1[-flds[[f]]] #train x
      X1_test<-X1[flds[[f]]] #train x
      X2_train<-X2[-flds[[f]]] #train x
      X2_test<-X2[flds[[f]]] #train x
      #pkyrs_train<-pkyrs[-flds[[f]]] #train x
      #pkyrs_test<-pkyrs[flds[[f]]] #train x
      Y_train=Y[-flds[[f]]] #train y
      Y_test=Y[flds[[f]]] #test y
      
      
      #kmeans clustering
      n.k<-NbClust(data=Gs2_train,min.nc = 2, max.nc = 9,method="kmeans",index="ch")$Best.nc["Number_clusters"]
      cl<-kmeans(Gs2_train,n.k)
      grp_assign_train<-cl$cluster #group assignement
      
      #distance matrix
      dist.matrix=matrix(,nrow=nrow(Gs2_test),ncol=n.k)
      for(k in 1: n.k){
        miu_k<-cl$centers[k,]
        dist.matrix[,k]<-apply(Gs2_test,1,function(x) sum((x-miu_k)^2))
      }
      
      #group assignment of test data
      grp_assign_test=apply(dist.matrix, 1, which.min)
      #grp_assign_prd[flds[[f]]]=grp_assign_test
      
      #calculate the y_est
      Y_prd_test=Y_test
      for(k in 1: n.k){
        y_k_train=Y_train[grp_assign_train==k]
        X1_k_train=X1_train[grp_assign_train==k]
        X2_k_train=X2_train[grp_assign_train==k]
        #pkyrs_k_train=pkyrs_train[grp_assign_train==k]
        fit_k<-lm(y_k_train ~ X1_k_train + X2_k_train)
        y_k_test=Y_test[grp_assign_test==k]
        X1_k_test=X1_test[grp_assign_test==k]
        X2_k_test=X2_test[grp_assign_test==k]
        #pkyrs_k_test=pkyrs_test[grp_assign_test==k]
        Y_prd_test[grp_assign_test==k]=fit_k$coefficients %*% rbind(1, X1_k_test, X2_k_test)
      }
      onecv.res<-list(pvalues=0,K=n.k, fld_id=flds[[f]], 
                      grp_assign=grp_assign_test, Y_prd_test=Y_prd_test)
      cvres.scluster[[f]]<-onecv.res
    }
    
    Y_prd<-Y
    grp_assign<-Y
    for(i in 1:length(flds)){
      id=cvres.scluster[[i]]$fld_id
      Y_prd[id]<-cvres.scluster[[i]]$Y_prd_test
      grp_assign[id]<-cvres.scluster[[i]]$grp_assign
    }
    
    RMSE[o]=sqrt(sum((Y_prd-Y)^2)/n)
  }
  
  final.ngs<-ngs[which.min(RMSE)]
  Gs2<-Gs[,order(ps)[1:final.ngs]]
  
  #kmeans clustering
  n.k<-NbClust(data=Gs2,min.nc = 2, max.nc = 9,method="kmeans",index="ch")$Best.nc["Number_clusters"]
  cl<-kmeans(Gs2,n.k)
  cluster<-cl$cluster
  return(list(n.k, order(ps)[1:final.ngs] ,cluster))
}

save(res.scluster, gen.res, file="scluster_K3_1000genes_9groups_exp0.5_gamma1_delta2.rdata")
library(mclust)
est.K<-c()
ARI<-c()
for(i in 1:100){
  est.K[i]<-res.scluster[[i]][[1]]
  ARI[i]<-adjustedRandIndex(res.scluster[[i]][[3]],gen.res$grp_assign[i,])
}
table(est.K)
ARI
median(ARI) #[1] 0.4166601

# Genes: FNs FPs
TPs<-c()
TNs<-c()
FNs<-c()
FPs<-c()
for(i in 1:100){
  TPs[i]<-sum(res.scluster[[i]][[2]] %in% 1:15)
  TNs[i]<-sum((1:1000)[-res.scluster[[i]][[2]]] %in% 16:1000)
  FNs[i]<-sum((1:1000)[-res.scluster[[i]][[2]]] %in% 1:15)
  FPs[i]<-sum(res.scluster[[i]][[2]] %in% 16:1000)
}

summary(FNs)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5       5       5       5       5       5 
summary(FPs)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0     0.0     2.5     5.5     5.0    20.0