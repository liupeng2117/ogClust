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
miu=c(1,6,11) # overall mean
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
    while(!is.na(abs(log_lik_up-log_lik)) & abs(log_lik_up-log_lik)>10^(-7) & iter <=max_iter){
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
    if(is.na(log_lik_up)) print("NA ll found!!")
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

res.pmbc <- foreach(rep = 1:n.rep, .packages = c("Brobdingnag")) %dopar% {
  print(paste("data", rep))
  df<-data[[rep]]
  G=as.matrix(df[,3:(NG+2)])
  G=apply(G,2,function(x) (x-mean(x))/sd(x))
  #Use BIC to select the number of clusters and tuning parameter lambda
  K.all=c(2,3,6,9)
  lambda.all=seq(6,14,2)
  BIC.all<-matrix(,ncol=length(lambda.all),nrow=length(K.all))
  for(l in 1:length(K.all)){
    K=K.all[l]
    #initialize using k means------
    data.kmeans<-kmeans(as.matrix(G),K)
    c_size<-data.kmeans$size
    pi_int<-c_size/n
    miu_int<-sapply(1:K,function(x) apply(as.matrix(G[data.kmeans$cluster==x,]),2,mean))
    dim(miu_int) #100 x K
    sigma_c<-sapply(1:K,function(x) apply(as.matrix(G[data.kmeans$cluster==x,]),2,var))
    sum_squares<-sapply(1:K,function(x) c_size[x]*sigma_c[,x])
    sigma_int<-apply(sum_squares,1,sum)/n
    
    #Fit mbc model and calculate BIC
    for(m in 1:length(lambda.all)){
      lambda=lambda.all[m]
      print(paste0("K=",K,", lambda=",lambda))
      BIC<-tryCatch({
        BIC<-em_mbc(t(G), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                    lambda=lambda, K=K, no_init=1, max_iter = 200)$BIC 
      }, error=function(e){
        BIC=Inf
        return(BIC)
      })
      BIC.all[l,m]<-BIC
    }
  }
  #select the best K and lambda based on BIC
  K=best.K=K.all[which(BIC.all==min(BIC.all,na.rm=T),arr.ind = TRUE)[1]] #select the best K
  best.lambda=lambda.all[which(BIC.all==min(BIC.all,na.rm=T),arr.ind = TRUE)[2]] #select the best lambda
  cl<-em_mbc(t(G), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                  lambda=best.lambda, K=K, no_init=1, max_iter = 200)
  theta_est<-cl$optimal_result
  mult_pdf<-matrix(,nrow=n,ncol=K)
  z_est<-matrix(,nrow=n, ncol=K)
  for(i in 1:n){
    for(j in 1:K){
      mult_pdf[i,j]<-as.numeric(mult_density(t(G)[,i],mu=theta_est$mu[,j],sigma=theta_est$sigma))
    }
  }
  
  for(i in 1:n){
    d<-brob(mult_pdf[i,])*theta_est$pi
    z_est[i,]<-as.numeric(d/sum(d))
  }
  
  return(list(best.K, best.lambda, theta_est, z_est))
}

save(res.pmbc,gen.res, file="pmbc_K3_1000genes_9groups_exp0.5_gamma1_delta5.rdata")

library(mclust)
est.K<-c()
ARI<-c()
for(i in 1:100){
  est.K[i]<-res.pmbc[[i]][[1]]
  ARI[i]<-adjustedRandIndex(apply(res.pmbc[[i]][[4]],1,which.max),gen.res$grp_assign[i,])
}
table(est.K)
ARI
median(ARI)

# Genes: FNs FPs
TPs<-c()
TNs<-c()
FNs<-c()
FPs<-c()
for(i in 1:100){
  TPs[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)>0) %in% 1:15)
  TNs[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)==0) %in% 16:1000)
  FNs[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)==0) %in% 1:15)
  FPs[i]<-sum(which(apply(abs(res.pmbc[[i]][[3]]$mu),1,mean)>0) %in% 16:1000)
}

summary(FNs)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#15      15      15      15      15      15 
summary(FPs)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#15      15      15      15      15      1


