rm(list=ls())
#source("~/subtyping/simulation/code/gen.data.K3.R")
library(survival)
library(reda)
source("survival_data_gen_peng_K3.R")
# Data generation
set.seed(123)
n.rep=100# repeat time 
NG=1000
n=600 # sample
np=2 # 2 covariates

beta1=0.5 # coefficient of clinical var1
beta2=-0.5 # coefficient of clinical var2
gamma1=c(3,3,3,3,3,0,0,0,0,0,-3,-3,-3,-3,-3) # 1st set coefficient of genes 1 to 15, others are 0
gamma2=c(0,0,0,0,0,3,3,3,3,3,-3,-3,-3,-3,-3) # 2nd set coefficient of genes 1 to 15, others are 0
miu=c(1,2,3) # overall mean
sigma=0.5 # overall variance
K=3
theta=list(beta1=beta1,
           beta2=beta2,
           offset=1,
           gamma1=gamma1,
           gamma2=gamma2,
           miu=miu,
           sigma=sigma)

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

gen.train<-survival_data_gen(n.rep,NG,n,np,theta,cluster_index)
data.train<-gen.train$data

gen.test<-survival_data_gen(n.rep,NG,n,np,theta,cluster_index)
data.test<-gen.test$data

standardize<-function(mat){
  mat.std<-apply(mat,2,function(x) (x-mean(x))/sd(x))
  return(mat.std)
}

#---------------------------------------------------#
#  Configuration Panel #
#--------------------------------------------------#
K.all=3 #Set K = truth
lambda.all<-seq(1.1,3,0.3) # penalty parameter
tune_pars=cbind(K=rep(K.all,each=length(lambda.all)),lambda=rep(lambda.all,length(K.all)))
#---------------------------------------------------#

library(doParallel) # foreach
cl<-makeCluster(10)
registerDoParallel(cl)

library(KLClust)
library(NbClust)
library(mclust)
res.skm <- foreach(rep = 1:n.rep, .packages=c("mclust","KLClust","NbClust","survival")) %dopar% {
  print(paste("data",rep))
  #train
  df.train<-data.train[[rep]]
  X1_train=df.train$X1
  X2_train=df.train$X2
  G_train=as.matrix(df.train[,3:(NG+2)])
  G_train=sapply(1:ncol(G_train), function(x) (G_train[,x]-mean(G_train[,x]))/sd(G_train[,x]))
  Y_train=df.train$time
  delta_train=df.train$event
  
  #test
  df.test<-data.test[[rep]]
  #   T_test=df.test$T
  X1_test=df.test$X1
  X2_test=df.test$X2
  G_test=as.matrix(df.test[,3:(NG+2)])
  G_test=sapply(1:ncol(G_test), function(x) (G_test[,x]-mean(G_test[,x]))/sd(G_test[,x]))
  Y_test=df.test$time
  delta_test=df.test$event
  
  onerep.res <- list()
  for (l in 1:nrow(tune_pars)) {
    K<-tune_pars[l, "K"]
    lambda <- tune_pars[l, "lambda"]
    cltrain<-KMeansSparseCluster(G_train,K=K,wbounds=lambda)
    #cltrain<-KMeansSparseCluster(Gs,K=kl$K,wbounds=kl$lambda)
    
    weight<-cltrain[[1]]$ws
    cluster<-cltrain[[1]]$Cs
    grp_assign_train<-cltrain[[1]]$Cs
    ARI <- adjustedRandIndex(grp_assign_train, gen.train$grp_assign[rep, ])
    #K=length(unique(grp_assign_train))
    n_genes_selected=sum(weight>0)
    n_true_genes_selected=sum(weight[1:15]>0)
    
    #distance matrix
    dist.matrix=matrix(,nrow=nrow(G_test),ncol=K)
    for(k in 1: K){
      miu_k<-apply(G_train[grp_assign_train==k,],2,mean)
      dist.matrix[,k]<-apply(G_test,1,function(x) weight %*% (x-miu_k)^2)
    }
    #group assignment of test data
    grp_assign_test=apply(dist.matrix, 1, which.min)
    ARI_test <- adjustedRandIndex(grp_assign_test, gen.test$grp_assign[rep, ])
    #calculate the y_est
    y_prd_test=Y_test
    for(k in 1: K){
      y_k_train=Y_train[grp_assign_train==k]
      delta_k_train=delta_train[grp_assign_train==k]
      x1_k_train=X1_train[grp_assign_train==k]
      x2_k_train=X2_train[grp_assign_train==k]
      #fit_k<-lm(y_k_train ~ x1_k_train + x2_k_train)
      fit_k<-survreg(Surv(y_k_train,delta_k_train) ~ x1_k_train + x2_k_train,dist = "loglogistic", robust=TRUE)
      y_k_test=Y_test[grp_assign_test==k]
      x1_k_test=X1_test[grp_assign_test==k]
      x2_k_test=X2_test[grp_assign_test==k]
      y_prd_test[grp_assign_test==k]=exp(fit_k$coefficients %*% rbind(1, x1_k_test, x2_k_test))
    }
    y_prd_test<-ifelse(y_prd_test>100,100,y_prd_test)
    
    #calculate RMSE
    RMSE=sqrt(sum((Y_test-y_prd_test)^2)/n)
    R2=1-sum((Y_test-y_prd_test)^2)/sum((Y_test-mean(Y_test))^2)
    
    onerep.res[[l]]<-list(K = K, best.lambda = lambda, 
                          grp_assign_train = grp_assign_train, grp_assign_test = grp_assign_test,
                          ARI = ARI, ARI_test = ARI_test, Y_prd_test=y_prd_test,
                          n_genes_selected=n_genes_selected, n_true_genes_selected=n_true_genes_selected)
    
  }
  return(onerep.res)
}

save(res.skm,gen.train,gen.test,file="survival_skmeans_K3_delta1_gamma3.rdata")

ARI <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
n_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
n_true_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2 <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)

for(j in 1:nrow(tune_pars)){
  for (i in 1:n.rep) {
    #train
    Y <- gen.train$data[[i]]$time
    Z <- gen.train$data[[i]]$Z
    X1 <- gen.train$data[[i]]$X1
    X2 <- gen.train$data[[i]]$X2
    #Y_prd <- res.skm[[i]][[j]]$Y_prd
    #e=Y-(miu[Z]+beta1*X1+beta2*X2)
    #RMSE[j,i]<- sqrt(sum((Y_prd - (Y-e))^2) / n)
    #R2[j,i] <- 1 - sum((Y_prd-(Y-e))^2) / sum(((Y-e) - mean(Y))^2)
    ARI[j,i] <- res.skm[[i]][[j]]$ARI
    n_markers[j,i]<- res.skm[[i]][[j]]$n_genes_selected
    n_true_markers[j,i]<- res.skm[[i]][[j]]$n_true_genes_selected    
    grp_assign_train <- res.skm[[i]][[j]]$grp_assign_train
    # pvalues[i]=kruskal.test(Y ~ grp_assign_train)$p.value
    print(table(grp_assign_train, gen.train$grp_assign[i,]))
    
    #test
    Y_test <- gen.test$data[[i]]$time
    Z_test <- gen.test$data[[i]]$Z
    X1_test <- gen.test$data[[i]]$X1
    X2_test <- gen.test$data[[i]]$X2
    Y_prd_test <- res.skm[[i]][[j]]$Y_prd_test
    RMSE_test[j,i]<- sqrt(sum((Y_prd_test - Y_test)^2) / n)
    R2_test[j,i] <- 1 - sum((Y_prd_test-Y_test)^2) / sum((Y_test - mean(Y_test))^2)
    ARI_test[j,i] <- res.skm[[i]][[j]]$ARI_test
  }
}

ARI 
#[1] 0.8113398 0.8612458 0.8187639 0.8350720 0.7644884
#[6] 0.8353246 0.8441081 0.4216053 0.7630068 0.4806123
mean(ARI) #0.7435567
# summary(pvalues)
ARI_test

mean(ARI_test)

n_markers

n_true_markers

#apply(RMSE,1,mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.072   2.188   2.265   2.409   2.312   3.235  
#apply(R2,1,mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3327  0.5072  0.5291  0.5090  0.5641  0.6020 

apply(RMSE_test,1,mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.072   2.188   2.265   2.409   2.312   3.235  
apply(R2_test,1,mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3327  0.5072  0.5291  0.5090  0.5641  0.6020 