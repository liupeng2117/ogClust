rm(list=ls())
source("~/subtyping/simulation/code/gen.data.K3.R")
#source("D:/research/Disease subtyping/simulation/gen.data.K3.R")
library(OGClust)
#library(mclust)
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

# Gaussian mixture model with subject specific probability --------------------------------------------------
standardize<-function(mat){
  mat.std<-apply(mat,2,function(x) (x-mean(x))/sd(x))
  return(mat.std)
}

#start<-proc.time()

library(doParallel) # foreach
cl<-makeCluster(10)
registerDoParallel(cl)

K.all=c(2,3,6,9) #Set K = truth
lambda.all<-c(0.02,0.03,0.05,0.1,0.15) # penalty parameter
tune_pars=cbind(K=rep(K.all,each=length(lambda.all)),lambda=rep(lambda.all,length(K.all)))

robust="none" # robust method
alpha=0.3 # for elastic net
# 5 fold
library(caret)
set.seed(123)
flds <- createFolds(1:n, k = 5, list = TRUE, returnTrain = FALSE)
#---------------------------------------------------#

#5 fold cross validation
cvres.all<-foreach(rep=1:n.rep,.packages = c("mclust","OGClust")) %dopar% {
  #print(rep)
  df<-data[[rep]]
  X1=df$X1
  X2=df$X2
  X=cbind(X1,X2)
  df[,3:(NG+2)]<-standardize(as.matrix(df[,3:(NG+2)])) #standardize
  G=as.matrix(df[,3:(NG+2)])
  Y=df$Y
  cvres.one<-list()
  
  for(f in 1:length(flds)){
    
    #train and test data
    df_train=df[-flds[[f]],] #train data
    df_test=df[flds[[f]],] #test data
    G_train=G[-flds[[f]],] #train g
    G_test=G[flds[[f]],] #test g
    X_train=X[-flds[[f]],] #train x
    X_test=X[flds[[f]],] #test x
    Y_train=Y[-flds[[f]]] #train y
    Y_test=Y[flds[[f]]] #test y
    
    res.all<-list()
    for(l in 1:nrow(tune_pars)){
      #parameters
      K=tune_pars[l,"K"]
      lambda=tune_pars[l,"lambda"]
      #2 initials and choose the best
      n.int=3
      est.par=matrix(ncol=((K-1)*(NG+1)+np+K+6),nrow=n.int)
      for(i in 1:n.int){
        beta_int <- runif(np, 0, 3)
        gamma_int <- runif((K - 1) * (NG + 1), 0, 1)
        beta0_int <- runif(K, 0, 3)
        sigma2_int <- runif(1, 1, 3)
        theta_int <- c(beta_int, gamma_int, beta0_int, sigma2_int)
        fit.res <- fit.OGClust(
          n = n, K = K, np = np, NG = NG, lambda = lambda, alpha = alpha,
          G = G_train, Y = Y_train, X = X_train, theta_int = theta_int, robust = robust
        )
        est.par[i, ] <- fit.res$res
      }
      #choose the best par estimates which has the maximum ll
      best.par<-est.par[which.min(est.par[,(K-1)*(NG+1)+K+np+5]),]
      res.all[[l]]<-best.par
    }
    
    BIC.all<-c()
    for (l in 1:nrow(tune_pars)){
      K=tune_pars[l,"K"]
      BIC.all[l]<-res.all[[l]][(K-1)*(NG+1)+np+K+5]
    }
    
    print("All K and lambda DONE")
    # best K
    K=best.K<-tune_pars[which.min(BIC.all),"K"]
    best.lambda=tune_pars[which.min(BIC.all),"lambda"]
    print(paste0("best K=", best.K))
    print(paste0("best lambda=", best.lambda))
    
    best.par<-res.all[[which.min(BIC.all)]]
    print(best.par[c(1:2,104:106)])
    best.lambda=best.par[(K-1)*(NG+1)+K+8]
    beta1_est=best.par[1]
    beta2_est=best.par[2]
    beta_est=c(beta1_est,beta2_est)
    gamma_est=best.par[3:((K-1)*(NG+1)+2)]
    miu_est=best.par[((K-1)*(NG+1)+3):((K-1)*(NG+1)+2+K)]
    sigma2_est=best.par[((K-1)*(NG+1)+K+3)]
    
    gamma_est_matrix=matrix(gamma_est,ncol=K-1,byrow=T)
    gamma_est_matrix=cbind(gamma_est_matrix,0)
    
    pai_est=sapply(1:K, function(k) 
      exp(cbind(1,G_train) %*% gamma_est_matrix[,k,drop=F])/(rowSums(exp(cbind(1,G_train) %*% gamma_est_matrix))))
    #calculate the expected value of Z
    f_est<-sapply(1:K, function(x) (1/sqrt(2*pi*sigma2_est))*exp(-(Y_train-miu_est[x]-X_train %*% beta_est)^2/(2*sigma2_est)))
    Y_prd_train=apply(sapply(1:K, function(x) pai_est[,x]*(miu_est[x]+X_train %*% beta_est)),1,sum)
    #calculate the expected value of Z
    w_est=sapply(1:K, function(k) 
      (pai_est[,k]*f_est[,k])/diag(pai_est %*% t(f_est)))
    grp_assign_train<-apply(w_est,1,which.max)
    
    #test group assignment and y prediction
    pai_est=sapply(1:K, function(k) 
      exp(cbind(1,G_test) %*% gamma_est_matrix[,k,drop=F])/(rowSums(exp(cbind(1,G_test) %*% gamma_est_matrix))))
    #pai_est=cbind(pai_est,1-rowSums(pai_est))
    #pai_prd[[f]]<-pai_est
    Y_prd_test=apply(sapply(1:K, function(x) pai_est[,x]*(miu_est[x]+X_test %*% beta_est)),1,sum)
    
    f_est<-sapply(1:K, function(x) (1/sqrt(2*pi*sigma2_est))*exp(-(Y_test-miu_est[x]-X_test %*% beta_est)^2/(2*sigma2_est)))
    #calculate the expected value of Z
    w_est=sapply(1:K, function(k) 
      (pai_est[,k]*f_est[,k])/diag(pai_est %*% t(f_est)))
    grp_assign_test<-apply(w_est,1,which.max)
    
    res<-list(K=best.K, lambda=best.lambda, fld_id=flds[[f]], grp_assign_train=grp_assign_train,
              grp_assign_test=grp_assign_test, Y_prd_test=Y_prd_test, BIC=BIC.all, theta_est=best.par)
    cvres.one[[f]]<-res
  }  
  return(cvres.one)
}

save(cvres.all, gen.res, file="ogclust_rmse_K3_1000genes_9groups_gamma3_delta3_exp1.rdata")

library(mclust)
RMSE<-c()
R2<-c()
#ARI<-c()
Y_prd<-rep(0,n)
grp_assign_prd<-Y_prd
for(rep in 1:n.rep){
  cvres.one<-cvres.all[[rep]]
  for(i in 1:length(flds)){
    id=cvres.one[[i]]$fld_id
    Y_prd[id]<-cvres.one[[i]]$Y_prd_test
  }
  Y<-gen.res$data[[rep]]$Y
  #ARI[rep]<-adjustedRandIndex(grp_assign_prd[-out.id],gen.res$grp_assign[rep, -out.id])
  RMSE[rep]=sqrt(sum((Y_prd-Y)^2)/n)
  R2[rep]=1-sum((Y-Y_prd)^2)/sum((Y-mean(Y))^2)
}
#summary(ARI)
summary(RMSE)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.686   1.702   1.740   1.731   1.753   1.773 
summary(R2)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.6077  0.6259  0.6341  0.6342  0.6384  0.6624