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

gen.res<-gen.data.K3(n.rep, NG, n, theta, cluster_index, exp=1)
data<-gen.res$data

library(doParallel) # foreach
cl<-makeCluster(10)
registerDoParallel(cl)
#=============== Model based clustering ================#
library(Brobdingnag)# for avoid the small number equals 0 in calculation
# calculate log likelihood for subject j and cluster k
mult_density<-function(x,mu,sigma){
  
  sum<-sum(dnorm(x,mean=mu,sd=sqrt(sigma),log=TRUE))
  return(sum) 
}

# EM algorithm
em_mbc<-function(data,mu=0,sigma=1,lambda,c_center=NULL,v_int=NULL,pi_int=NULL,K=2,no_init=1,max_iter=200){
  result_list<-list()
  log_lik_vector<-c()
  unpen_lik_vector<-c()
  time<-c()
  #Column refers to samples and row refers to genes 
  p<-dim(data)[1] # number of variables
  n<-dim(data)[2] # number of subjects
  
  for(init in 1:no_init){
    start<-Sys.time()
    #========== E-step: initialization===========#
    
    if(length(c_center)==0){
      # assume the means for each cluster follows the distribution N(0,1)
      mu_int<-matrix(rnorm(p*K,mu,sigma),nrow=p,ncol=K)
      #mu[which(mu==0)]<-100
    } else{
      mu_int=c_center
    }
    
    v_int<-ifelse(length(v_int)!=p,rep(sigma,p),v_int) # initial value of variance of each variables
    pi_int<-ifelse(length(pi_int)!=K,rep(1/K,K),pi_int) # initial value of pi
    z_int<-matrix(,nrow=n,ncol=K) # initial value of prob in cluster k for each subject
    
    mult_pdf<-matrix(,nrow=n,ncol=K)
    for(i in 1:n){
      for(j in 1:K){
        mult_pdf[i,j]<-as.numeric(mult_density(data[,i],mu=mu_int[,j],sigma=v_int))
      }
    }
    
    for(i in 1:n){
      d<-brob(mult_pdf[i,])*pi_int
      z_int[i,]<-as.numeric(d/sum(d))
    }
    
    pi<-pi_int
    v<-v_int
    mu<-mu_int
    z<-z_int
    
    #========= M step ==========#
    iter=1
    log_lik<-1 # initialize
    log_lik_up<-0 #initialize
    while(abs(log_lik_up-log_lik)>10^(-7) & iter <=max_iter){
      # log likelihood before updating
      log_lik<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf[x,])))-lambda*sum(abs(mu))
      #update pi
      pi_up<-apply(z,2,sum)/n
      #pi_up<-ifelse(pi_up==0,1*10^(-16),pi_up)
      #update v
      update_v<-function(x){
        sig<-sum(sapply(1:K,function(y) sum(z[,y]*(data[x,]-mu[x,y])^2)))/n
        return(sig)
      }
      v_up<-sapply(1:p,function(x) update_v(x))
      
      #update mu
      mu_tu<-matrix(,nrow=p,ncol=K)
      mu_up<-matrix(,nrow=p,ncol=K)
      
      for(i in 1:K){
        temp<-sapply(1:n,function(x) z[x,i]*data[,x]) 
        mu_tu[,i]<-apply(temp,1,sum)/sum(z[,i])
        mu_up[,i]<-ifelse(lambda<=abs(apply(temp,1,sum)/v),sign(mu_tu[,i])*(abs(mu_tu[,i])-(lambda/sum(z[,i]))*v),0)
      }
      
      
      #update z
      z_up<-matrix(,nrow=n,ncol=K)
      
      for(i in 1:n){
        for(j in 1:K){
          mult_pdf[i,j]<-as.numeric(mult_density(data[,i],mu=mu_up[,j],sigma=v_up))
        }
      }
      
      for(i in 1:n){
        d<-brob(mult_pdf[i,])*pi_up
        z_up[i,]<-as.numeric(d/sum(d))
      }
      
      # update parameters values
      pi<-pi_up
      v<-v_up
      mu<-mu_up
      z<-z_up
      #log likelihood after updating
      log_lik_up<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf[x,])))-lambda*sum(abs(mu))
      unpen_lik_up<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf[x,])))
      iter=iter+1
    }
    
    end<-Sys.time()
    time[init]<-end-start
    result_list[[init]] <- list('mu'=mu,'sigma'=v,'log_lik'=log_lik,'z'=z,'pi'=pi,'no_init'=init) 
    log_lik_vector[init]<-log_lik_up 
    unpen_lik_vector[init]<-unpen_lik_up
    #print(init)
  }
  max_lik<-unpen_lik_vector[which.max(log_lik_vector)]
  
  # the optimal initials is the one with the maximum log likelihood
  optimal_result<-result_list[[which.max(log_lik_vector)]]  
  max_mu<-optimal_result$mu
  s<-apply(max_mu,1,function(x) sum(x==0))
  #d<-sum(s[which(s!=1)])
  q<-sum(s)
  #BIC
  if(lambda!=0) {
    BIC<--2*max_lik+log(n)*(K+p+K*p-1-q)
  } else if(lambda==0) {
    BIC<--2*max_lik+log(n)*(K+p+K*p-1)
  }
  
  list('optimal_result'=optimal_result,'result_list'=result_list,'time'=time,'BIC'=BIC)  
}



# penalized model based clustering
# cluster -----------------------------------------------------------------
library(caret)
flds <- createFolds(1:n, k = 5, list = TRUE, returnTrain = FALSE)

library(doParallel) # foreach
cl<-makeCluster(20)
registerDoParallel(cl)

#10-fold cross validation
res.pmbc<-foreach(rep=1:n.rep,.packages="Brobdingnag") %dopar% {
  print(rep)
  Y_prd=rep(0,n)
  grp_assign_prd=rep(0,n)
  #y_prd_cv_list<-foreach(f=1:10, .packages="Brobdingnag") %dopar% {
  y_prd_cv_list=list()
  for(f in 1:5){
    #train
    df.train<-data[[rep]][-flds[[f]],]
    #   T_train=df.train$T
    X1_train=df.train$X1
    X2_train=df.train$X2
    G_train=as.matrix(df.train[,3:(NG+2)])
    # standardized the data, such that each attribute has mean 0 and sd 1
    G_train<-apply(G_train,2,function(x) (x-mean(x))/sd(x))
    Y_train=df.train$Y
    
    #test
    df.test<-data[[rep]][flds[[f]],]
    #   T_test=df.test$T
    X1_test=df.test$X1
    X2_test=df.test$X2
    G_test=as.matrix(df.test[,3:(NG+2)])
    G_test<-apply(G_test,2,function(x) (x-mean(x))/sd(x))
    Y_test=df.test$Y
    
    #Use BIC to select the number of clusters and tuning parameter lambda
    K.all=c(2,3,6,9)
    lambda.all=seq(16,24,2)
    BIC.all<-matrix(,ncol=length(lambda.all),nrow=length(K.all))
    for(l in 1:length(K.all)){
      K=K.all[l]
      #initialize using k means------
      data.kmeans<-kmeans(as.matrix(G_train),K)
      c_size<-data.kmeans$size
      pi_int<-c_size/n
      miu_int<-sapply(1:K,function(x) apply(as.matrix(G_train[data.kmeans$cluster==x,]),2,mean))
      dim(miu_int) #100 x K
      sigma_c<-sapply(1:K,function(x) apply(as.matrix(G_train[data.kmeans$cluster==x,]),2,var))
      sum_squares<-sapply(1:K,function(x) c_size[x]*sigma_c[,x])
      sigma_int<-apply(sum_squares,1,sum)/n
      
      #Fit mbc model and calculate BIC
      for(m in 1:length(lambda.all)){
        lambda=lambda.all[m]
        print(paste0("K=",K,", lambda=",lambda))
        BIC<-tryCatch({
          BIC<-em_mbc(t(G_train), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                      lambda=lambda, K=K, no_init=1, max_iter = 200)$BIC 
        }, error=function(e){
          BIC=Inf
          return(BIC)
        })
        BIC.all[l,m]<-BIC
      }
    }
    #select the best K and lambda based on BIC
    K=best.K=K.all[which(BIC.all==min(BIC.all),arr.ind = TRUE)[1]] #select the best K
    best.lambda=lambda.all[which(BIC.all==min(BIC.all),arr.ind = TRUE)[2]] #select the best lambda
    cltrain<-em_mbc(t(G_train), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                    lambda=best.lambda, K=K, no_init=1, max_iter = 200)
    pi_est_train<-cltrain$optimal_result$pi
    mu_est_train<-cltrain$optimal_result$mu
    sigma_est_train<-cltrain$optimal_result$sigma
    grp_assign_train<-apply(cltrain$optimal_result$z,1,which.max)
    
    mult_pdf<-matrix(nrow=n,ncol=K)
    postprob.matrix<-matrix(,nrow=n, ncol=K)
    for(i in 1:n){
      for(j in 1:K){
        mult_pdf[i,j]<-as.numeric(mult_density(G_test[,i],mu=mu_est_train[,j],sigma=sigma_est_train))
      }
    }
    
    for(i in 1:n){
      d<-brob(mult_pdf[i,])*pi_est_train
      postprob.matrix[i,]<-as.numeric(d/sum(d))
    }
    #group assignment of test data
    grp_assign_test=apply(postprob.matrix, 1, which.max)
    grp_assign_prd[flds[[f]]]=grp_assign_test
    #calculate the y_est
    y_prd_cv=Y_test
    for(k in unique(grp_assign_test)){
      y_k_train=Y_train[grp_assign_train==k]
      x1_k_train=X1_train[grp_assign_train==k]
      x2_k_train=X2_train[grp_assign_train==k]
      fit_k<-lm(y_k_train ~ x1_k_train + x2_k_train)
      y_k_test=Y_test[grp_assign_test==k]
      x1_k_test=X1_test[grp_assign_test==k]
      x2_k_test=X2_test[grp_assign_test==k]
      y_prd_cv[grp_assign_test==k]=fit_k$coefficients %*% rbind(1, x1_k_test, x2_k_test)
    }
    onecv.res<-list(y_prd_cv, grp_assign_test)
    y_prd_cv_list[[f]]<-onecv.res
  }
  for(f in 1:5){
    Y_prd[flds[[f]]]<-y_prd_cv_list[[f]][[1]]
    grp_assign_prd[flds[[f]]]<-y_prd_cv_list[[f]][[2]]
  }
  #calculate RMSE
  RMSE=sqrt(sum((data[[rep]]$Y-Y_prd)^2)/n)
  R2=1-sum((data[[rep]]$Y-Y_prd)^2)/sum((data[[rep]]$Y-mean(data[[rep]]$Y))^2)
  return(list(RMSE=RMSE, R2=R2, grp_assign_prd=grp_assign_prd))
}

save(res.pmbc, gen.res, file="pmbc_K3_1000genes_9groups_exp1_gamma3_delta5.rdata")

library(mclust)
RMSE.all<-c()
R2<-c()
ARI.all<-c()
for(rep in 1:n.rep){
  RMSE.all[rep]<-res.pmbc[[rep]]$RMSE
  R2[rep]<-res.pmbc[[rep]]$R2
  grp_assign_prd<-res.pmbc[[rep]]$grp_assign_prd
  ARI.all[rep]<-adjustedRandIndex(gen.res$grp_assign[rep,], grp_assign_prd)
  print(table(gen.res$grp_assign[rep,], grp_assign_prd))
}
summary(RMSE.all)
summary(R2)
summary(ARI.all)
