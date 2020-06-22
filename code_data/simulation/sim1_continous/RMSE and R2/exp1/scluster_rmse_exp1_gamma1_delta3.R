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

gen.res<-gen.data.K3(n.rep, NG, n, theta, cluster_index, exp=1)
data<-gen.res$data

library(doParallel) # foreach
cl<-makeCluster(10)
registerDoParallel(cl)

# Sparse kmeans
# cluster -----------------------------------------------------------------
#library(sparcl)
library(caret)
#library(KLClust)
flds <- createFolds(1:n, k = 5, list = TRUE, returnTrain = FALSE)

#10-fold cross validation
res.scluster <- foreach(rep = 1:n.rep, .packages=c("NbClust","caret")) %dopar% {
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
    #kl<-KL.S4(t(G_train),lambda_list = list(seq(4,15,1),seq(4,14,1),seq(4,15,1)), 
    #      k_vector=c(2,3,6),num.cores = 1)
    
    # select the outcome related genes
    ps<-c()
    for(i in 1:ncol(G_train)){
      ps[i]<-summary(lm(Y_train ~ G_train[,i]))$coefficients[2,4]
    }
    ngs<-c(10,15,30,60,100)
    RMSE=c()
    # cross validation to select the best threshold
    
    for(o in 1:length(ngs)){
      ngsm<-ngs[o]
      G_train2<-G_train[,order(ps)[1:ngsm]]
      
      # 5 fold cross validation
      #parallel computing
      #library(doParallel) # foreach
      #cl<-makeCluster(4)
      #registerDoParallel(cl)
      
      # 5 fold
      n1=nrow(G_train2)
      set.seed(123)
      flds2 <- createFolds(1:n1, k = 5, list = TRUE, returnTrain = FALSE)
      cvres.scluster<-list()
      for(f2 in 1:length(flds2)){
        #train and test data
        G_train3<-G_train2[-flds2[[f2]],] #train g
        G_test3<-G_train2[flds2[[f2]],] #test g
        X1_train3<-X1_train[-flds2[[f2]]] #train x
        X1_test3<-X1_train[flds2[[f2]]] #train x
        X2_train3<-X2_train[-flds2[[f2]]] #train x
        X2_test3<-X2_train[flds2[[f2]]] #train x
        Y_train3=Y_train[-flds2[[f2]]] #train y
        Y_test3=Y_train[flds2[[f2]]] #test y
        
        
        #kmeans clustering
        n.k<-NbClust(data=G_train3,min.nc = 2, max.nc = 9,method="kmeans",index="ch")$Best.nc["Number_clusters"]
        cl<-kmeans(G_train3,n.k)
        grp_assign_train3<-cl$cluster #group assignement
        
        #distance matrix
        dist.matrix=matrix(,nrow=nrow(G_test3),ncol=n.k)
        for(k in 1: n.k){
          miu_k<-cl$centers[k,]
          dist.matrix[,k]<-apply(G_test3,1,function(x) sum((x-miu_k)^2))
        }
        
        #group assignment of test data
        grp_assign_test3=apply(dist.matrix, 1, which.min)
        #grp_assign_prd[flds2[[f]]]=grp_assign_test3
        
        #calculate the y_est
        Y_prd_test3=Y_test3
        for(k in 1: n.k){
          y_k_train3=Y_train3[grp_assign_train3==k]
          X1_k_train3=X1_train3[grp_assign_train3==k]
          X2_k_train3=X2_train3[grp_assign_train3==k]
          #pkyrs_k_train3=pkyrs_train3[grp_assign_train3==k]
          fit_k<-lm(y_k_train3 ~ X1_k_train3 + X2_k_train3)
          y_k_test3=Y_test3[grp_assign_test3==k]
          X1_k_test3=X1_test3[grp_assign_test3==k]
          X2_k_test3=X2_test3[grp_assign_test3==k]
          #pkyrs_k_test3=pkyrs_test3[grp_assign_test3==k]
          Y_prd_test3[grp_assign_test3==k]=fit_k$coefficients %*% rbind(1, X1_k_test3, X2_k_test3)
        }
        onecv.res<-list(pvalues=0,K=n.k, fld_id=flds2[[f2]], 
                        grp_assign=grp_assign_test3, Y_prd_test3=Y_prd_test3)
        cvres.scluster[[f2]]<-onecv.res
      }
      
      Y_prd2<-data[[rep]]$Y
      #grp_assign<-Y
      for(i in 1:length(flds2)){
        id=cvres.scluster[[i]]$fld_id
        Y_prd2[id]<-cvres.scluster[[i]]$Y_prd_test
       # grp_assign[id]<-cvres.scluster[[i]]$grp_assign
      }
      
      RMSE[o]=sqrt(sum((Y_prd2-data[[rep]]$Y)^2)/n)
    }
    
    final.ngs<-ngs[which.min(RMSE)]
    G_train2<-G_train[,order(ps)[1:final.ngs]]
    G_test2<-G_test[,order(ps)[1:final.ngs]]
    
    #kmeans clustering
    n.k<-NbClust(data=G_train2,min.nc = 2, max.nc = 6,method="kmeans",index="ch")$Best.nc["Number_clusters"]
    cltrain<-kmeans(G_train2,n.k)
    grp_assign_train<-cltrain$cluster
    
    #distance matrix
    dist.matrix=matrix(,nrow=nrow(G_test2),ncol=n.k)
    for(k in 1: n.k){
      miu_k<-cltrain$centers[k,]
      dist.matrix[,k]<-apply(G_test2,1,function(x) sum((x-miu_k)^2))
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

save(res.scluster, gen.res, file="scluster_rmse_K3_1000genes_9groups_exp1_gamma1_delta3.rdata")

library(mclust)
RMSE.all<-c()
R2<-c()
ARI.all<-c()
for(rep in 1:n.rep){
  RMSE.all[rep]<-res.scluster[[rep]]$RMSE
  R2[rep]<-res.scluster[[rep]]$R2
  grp_assign_prd<-res.scluster[[rep]]$grp_assign_prd
  ARI.all[rep]<-adjustedRandIndex(gen.res$grp_assign[rep,], grp_assign_prd)
  print(table(gen.res$grp_assign[rep,], grp_assign_prd))
}
summary(RMSE.all)
summary(R2)
summary(ARI.all)
