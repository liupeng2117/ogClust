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



##lung disease##
load("/home/pel67/subtyping/lung/data/LungData.all.Rdata")
#library(KLClust)
gdata<-LungData.all$gene
dim(gdata) # 319 x 15966

#filter out genes with low mean count and low sd
g.mean<-apply(gdata,2,mean)
cut.mean=quantile(g.mean,probs=0.25)
gdata=gdata[,g.mean>cut.mean] # remove genes with 25% lowest variance or 25% lowest mean expression  

g.sd=apply(gdata,2, sd) #cv in original scale
cut.sd=sort(g.sd,decreasing = T)[2000]
gdata=gdata[,g.sd>=cut.sd]
dim(gdata) #319 x 2000

#standardize gene expression
gdata=sapply(1:ncol(gdata), function(x) (gdata[,x]-mean(gdata[,x]))/sd(gdata[,x]))

# exlude patients with missing outcome
# exlude genes if correlation with outcome is close to 0(<=0.2)
index.Y=!is.na(LungData.all$clinic[,"fev1pd1a"]) #exclude if outcome is missing
Y=LungData.all$clinic[index.Y,"fev1pd1a"] 
#r<-c()
#for(i in 1:ncol(gdata)){
#  r[i]<-cor(Y,gdata[index.Y,i])
#}
#index.G<-abs(r)>0.4
gdata2=gdata[index.Y,]
dim(gdata2) #234 x 2000

#covariates to adjust for 
X1=age=LungData.all$clinic[index.Y,"AGE"] 
X2=bmi=LungData.all$clinic[index.Y,"BMI"] 
X3=gender=LungData.all$clinic[index.Y,"GENDER"] #gender

# 5 fold cross validation
#parallel computing
library(doParallel) # foreach
cl<-makeCluster(5)
registerDoParallel(cl)

# 5 fold
n=nrow(gdata2)
library(caret)
set.seed(123)
flds <- createFolds(1:n, k = 5, list = TRUE, returnTrain = FALSE)

cvres.skm<-foreach(f=1:length(flds),.packages = "Brobdingnag") %dopar% {
  #train and test data
  gdata2_train<-gdata2[-flds[[f]],] #train g
  gdata2_test<-gdata2[flds[[f]],] #test g
  age_train<-age[-flds[[f]]] #train x
  age_test<-age[flds[[f]]] #train x
  bmi_train<-bmi[-flds[[f]]] #train x
  bmi_test<-bmi[flds[[f]]] #train x
  gender_train<-gender[-flds[[f]]] #train x
  gender_test<-gender[flds[[f]]] #train x
  Y_train=Y[-flds[[f]]] #train y
  Y_test=Y[flds[[f]]] #test y

  
  #------ choose the optimal k and lambda ---------#
  K=3
  lambda.all=seq(16,24,2)
  BIC.all<-matrix(,ncol=length(lambda.all),nrow=1)
  
  #initialize using k means------
  data.kmeans<-kmeans(as.matrix(gdata2_train),K)
  c_size<-data.kmeans$size
  pi_int<-c_size/n
  miu_int<-sapply(1:K,function(x) apply(as.matrix(gdata2_train[data.kmeans$cluster==x,]),2,mean))
  dim(miu_int) #100 x K
  sigma_c<-sapply(1:K,function(x) apply(as.matrix(gdata2_train[data.kmeans$cluster==x,]),2,var))
  sum_squares<-sapply(1:K,function(x) c_size[x]*sigma_c[,x])
  sigma_int<-apply(sum_squares,1,sum)/n
  
  #Fit mbc model and calculate BIC
  for(m in 1:length(lambda.all)){
    lambda=lambda.all[m]
    print(paste0("K=",K,", lambda=",lambda))
    BIC<-tryCatch({
      BIC<-em_mbc(t(gdata2_train), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
                  lambda=lambda, K=K, no_init=1, max_iter = 200)$BIC 
    }, error=function(e){
      BIC=Inf
      return(BIC)
    })
    BIC.all[1,m]<-BIC
  }
  #select the best K and lambda based on BIC
  #K=best.K=K.all[which(BIC.all==min(BIC.all,na.rm=T),arr.ind = TRUE)[1]] #select the best K
  K=best.K=3
  best.lambda=lambda.all[which(BIC.all==min(BIC.all,na.rm=T),arr.ind = TRUE)[2]] #select the best lambda
  
  #pmbc clustering
  cl_train<-em_mbc(t(gdata2_train), c_center=miu_int, v_int=sigma_int, pi_int=pi_int,
             lambda=best.lambda, K=best.K, no_init=1, max_iter = 200)
  postprob<-cl_train$optimal_result$z # weight
  grp_assign_train<-apply(postprob,1,which.max) #group assignement
  
  #test the association between cluster assignment and clinical variables
  pvalues<-c()
  for(i in 1:ncol(LungData.all$clinic)){
    values=LungData.all$clinic[index.Y,i][-flds[[f]]]
    tempdata=cbind(clusters=grp_assign_train,values=values)
    pvalues[i]<-kruskal.test(values ~ clusters, data = tempdata)$p.value 
  }
  names(pvalues)<-colnames(LungData.all$clinic)
  
  #distance matrix
  weight=rep(1,ncol(gdata2_test))
  dist.matrix=matrix(,nrow=nrow(gdata2_test),ncol=best.K)
  for(k in 1: best.K){
    miu_k<-apply(gdata2_train[grp_assign_train==k,],2,mean)
    dist.matrix[,k]<-apply(gdata2_test,1,function(x) weight %*% (x-miu_k)^2)
  }
  
  #group assignment of test data
  grp_assign_test=apply(dist.matrix, 1, which.min)
  #grp_assign_prd[flds[[f]]]=grp_assign_test
  
  #calculate the y_est
  Y_prd_test=Y_test
  for(k in 1: best.K){
    y_k_train=Y_train[grp_assign_train==k]
    age_k_train=age_train[grp_assign_train==k]
    bmi_k_train=bmi_train[grp_assign_train==k]
    gender_k_train=gender_train[grp_assign_train==k]
    fit_k<-lm(y_k_train ~ age_k_train + bmi_k_train + gender_k_train)
    y_k_test=Y_test[grp_assign_test==k]
    age_k_test=age_test[grp_assign_test==k]
    bmi_k_test=bmi_test[grp_assign_test==k]
    gender_k_test=gender_test[grp_assign_test==k]
    Y_prd_test[grp_assign_test==k]=fit_k$coefficients %*% rbind(1, age_k_test, bmi_k_test, gender_k_test)
  }
  onecv.res<-list(pvalues=pvalues,K=best.K,lambda=best.lambda, fld_id=flds[[f]], 
                  grp_assign=grp_assign_test, Y_prd_test=Y_prd_test)
  return(onecv.res)
}
save.image(file="lung_disease_pmbc_cv_top2000_K3.Rdata")

Y_prd<-Y
grp_assign<-Y
for(i in 1:length(flds)){
  id=cvres.skm[[i]]$fld_id
  Y_prd[id]<-cvres.skm[[i]]$Y_prd_test
  grp_assign[id]<-cvres.skm[[i]]$grp_assign
}

RMSE=sqrt(sum((Y_prd-Y)^2)/n)
RMSE
R2=1-sum((Y-Y_prd)^2)/sum((Y-mean(Y))^2)
R2

grp_diag<-ifelse(LungData.all$diagnosis$diagnosis=="COPD",1,2)[index.Y]
table(grp_diag, grp_assign)

pvalues=kruskal.test(Y ~ grp_assign)$p.value
pvalues
#1.083167e-06