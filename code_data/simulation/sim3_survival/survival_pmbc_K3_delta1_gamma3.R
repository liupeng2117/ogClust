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
lambda.all<-seq(16,24,2) # penalty parameter
tune_pars=cbind(K=rep(K.all,each=length(lambda.all)),lambda=rep(lambda.all,length(K.all)))
#---------------------------------------------------#

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


library(doParallel) # foreach
cl<-makeCluster(20)
registerDoParallel(cl)

#10-fold cross validation
res.pmbc<-foreach(rep=1:n.rep,.packages=c("Brobdingnag","survival","mclust")) %dopar% {
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
    
    #initialize using k means------
    data.kmeans<-kmeans(as.matrix(G_train),K)
    c_size<-data.kmeans$size
    pi_int<-c_size/n
    miu_int<-sapply(1:K,function(x) apply(as.matrix(G_train[data.kmeans$cluster==x,]),2,mean))
    dim(miu_int) #100 x K
    sigma_c<-sapply(1:K,function(x) apply(as.matrix(G_train[data.kmeans$cluster==x,]),2,var))
    sum_squares<-sapply(1:K,function(x) c_size[x]*sigma_c[,x])
    sigma_int<-apply(sum_squares,1,sum)/n
    
    cltrain<-em_mbc(t(G_train), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                    lambda=lambda, K=K, no_init=1, max_iter = 200)
    pi_est_train<-cltrain$optimal_result$pi
    mu_est_train<-cltrain$optimal_result$mu
    sigma_est_train<-cltrain$optimal_result$sigma
    grp_assign_train<-apply(cltrain$optimal_result$z,1,which.max)
    n_genes_selected=sum(apply(mu_est_train,1,function(x) sum(abs(x)))>0)
    n_true_genes_selected=sum(apply(mu_est_train,1,function(x) sum(abs(x)))[1:15]>0)
    ARI <- adjustedRandIndex(grp_assign_train, gen.train$grp_assign[rep, ])
    mult_pdf<-matrix(nrow=n,ncol=K)
    postprob.matrix<-matrix(nrow=n, ncol=K)
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

save(res.pmbc, gen.train, gen.test, file="survival_pmbc_K3_delta1_gamma3.rdata")

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
    #Y_prd <- res.pmbc[[i]][[j]]$Y_prd
    #e=Y-(miu[Z]+beta1*X1+beta2*X2)
    #RMSE[j,i]<- sqrt(sum((Y_prd - (Y-e))^2) / n)
    #R2[j,i] <- 1 - sum((Y_prd-(Y-e))^2) / sum(((Y-e) - mean(Y))^2)
    ARI[j,i] <- res.pmbc[[i]][[j]]$ARI
    n_markers[j,i]<- res.pmbc[[i]][[j]]$n_genes_selected
    n_true_markers[j,i]<- res.pmbc[[i]][[j]]$n_true_genes_selected    
    grp_assign_train <- res.pmbc[[i]][[j]]$grp_assign_train
    # pvalues[i]=kruskal.test(Y ~ grp_assign_train)$p.value
    print(table(grp_assign_train, gen.train$grp_assign[i,]))
    
    #test
    Y_test <- gen.test$data[[i]]$time
    Z_test <- gen.test$data[[i]]$Z
    X1_test <- gen.test$data[[i]]$X1
    X2_test <- gen.test$data[[i]]$X2
    Y_prd_test <- res.pmbc[[i]][[j]]$Y_prd_test
    RMSE_test[j,i]<- sqrt(sum((Y_prd_test - Y_test)^2) / n)
    R2_test[j,i] <- 1 - sum((Y_prd_test-Y_test)^2) / sum((Y_test - mean(Y_test))^2)
    ARI_test[j,i] <- res.pmbc[[i]][[j]]$ARI_test
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