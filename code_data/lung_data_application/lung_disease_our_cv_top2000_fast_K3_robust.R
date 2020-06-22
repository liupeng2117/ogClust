##lung disease##
library(OGClust)
#source("/home/pel67/subtyping/lung/code/EM_function_fast_elasticnet.R")
#source("/home/pel67/subtyping/lung/code/fit.mixture.model_function_elasticnet.R")
load("/home/pel67/subtyping/lung/data/LungData.all.Rdata")

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
index.Y=!is.na(LungData.all$clinic[,"fev1pd1a"]) #exclude if outcome is missing
Y=LungData.all$clinic[index.Y,"fev1pd1a"]
gdata2=gdata[index.Y,]
dim(gdata2) #234 x 2000

#sapply(1:50, function(x) cor(Y,gdata2[,x]))

#covariates to adjust for 
X1=age=LungData.all$clinic[index.Y,"AGE"] 
X2=bmi=LungData.all$clinic[index.Y,"BMI"] 
X3=gender=LungData.all$clinic[index.Y,"GENDER"] #gender
G=gdata2
df=data.frame(X1=age,X2=bmi,X3=gender, G,Y=Y)

# 5 fold cross validation
#parallel computing
library(doParallel) # foreach
cl<-makeCluster(21)
registerDoParallel(cl)

K.all=3
#lambda.all=c(2,4,6,8,10)
lambda.all=c(0.13, 0.11, 0.097834296, 0.051010883, 0.026597117)
#tune_pars=cbind(K=rep(K.all,each=length(lambda.all)),lambda=rep(lambda.all,length(K.all)))

n=nrow(gdata2)

# 5 fold
library(caret)
set.seed(123)
flds <- createFolds(1:n, k = 5, list = TRUE, returnTrain = FALSE)

#cross validatoin
cvres.our1000.all<-list()
for(ii in 1:length(lambda.all)){
  tune_pars=cbind(K=rep(K.all,each=length(lambda.all[ii])),lambda=rep(lambda.all[ii],length(K.all)))
  cvres.our1000<-list()
  for(f in 1:length(flds)){
    print(paste0("fold",f))
    #train and test data
    df_train=df[-flds[[f]],] #train data
    df_test=df[flds[[f]],] #test data
    G_train<-G[-flds[[f]],] #train g
    G_test<-G[flds[[f]],] #test g
    age_train<-age[-flds[[f]]] #train age
    age_test<-age[flds[[f]]] #train age
    bmi_train<-bmi[-flds[[f]]] #train BMI
    bmi_test<-bmi[flds[[f]]] #train BMI
    gender_train<-gender[-flds[[f]]] #train package years
    gender_test<-gender[flds[[f]]] #train package years
    X_train=cbind(age_train,bmi_train,gender_train) #train x
    X_test=cbind(age_test,bmi_test,gender_test) #test x
    Y_train=Y[-flds[[f]]] #train y
    Y_test=Y[flds[[f]]] #test y
    
    # number of genes and patients
    NG=ncol(G_train)
    n=nrow(G_train)
    np=ncol(X_train)
    
    res.all<-foreach(l=1:nrow(tune_pars), .packages = "OGClust") %dopar% {
      #parameters
      K=tune_pars[l,"K"]
      lambda=tune_pars[l,"lambda"]
      res.oneK=matrix(,nrow=length(lambda.all),ncol=((K-1)*(NG+1)+np+K+6))
      #2 initials and choose the best
      n.int=2
      est.par=matrix(,ncol=((K-1)*(NG+1)+np+K+6),nrow=n.int)
      #theta.int.all=matrix(,ncol=(K-1)*(NG+1)+np+K+1,nrow=n.int)
      for(i in 1:n.int){
        beta_int <- runif(np, 0, 3)
        gamma_int <- runif((K - 1) * (NG + 1), 0, 1)
        beta0_int <- runif(K, 0, 3)
        sigma2_int <- runif(1, 1, 3)
        theta_int <- c(beta_int, gamma_int, beta0_int, sigma2_int)
        fit.res<-fit.OGClust(n=n, K=K, np=np, NG=NG, lambda=lambda, 
                             alpha=0.5, G=G, Y=Y, X=X, theta_int = theta_int,
                             robust="hubertf2")
        est.par[i,]<-fit.res$res
      }
      #choose the best par estimates which has the maximum ll
      best.par<-est.par[which.min(est.par[,(K-1)*(NG+1)+K+np+5]),]
      return(best.par)
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
    print(res.all[[which.min(BIC.all)]])
    best.par=res.all[[which.min(BIC.all)]]
    
    #choose the best par estimates 
    beta_est=best.par[1:(np)]
    gamma_est=best.par[(np+1):((K-1)*(NG+1)+np)]
    beta0_est=best.par[((K-1)*(NG+1)+np+1):((K-1)*(NG+1)+K+np)]
    sigma2_est=best.par[(K-1)*(NG+1)+np+K+1]
    gamma_est_matrix=matrix(gamma_est,ncol=K-1,byrow=TRUE)
    gamma_est_matrix=cbind(gamma_est_matrix,0)
    
    #train group assignment and y prediction
    pai_est=sapply(1:K, function(k) 
      exp(cbind(1,G_train) %*% gamma_est_matrix[,k,drop=F])/(rowSums(exp(cbind(1,G_train) %*% gamma_est_matrix))))
    f_est<-sapply(1:K, function(x) (1/sqrt(2*pi*sigma2_est))*exp(-(Y_train-beta0_est[x]-X_train %*% beta_est)^2/(2*sigma2_est)))
    Y_prd_train=apply(sapply(1:K, function(x) pai_est[,x]*(beta0_est[x]+X_train %*% beta_est)),1,sum)
    #calculate the expected value of Z
    w_est=sapply(1:K, function(k) 
      (pai_est[,k]*f_est[,k])/diag(pai_est %*% t(f_est)))
    grp_assign_train<-apply(w_est,1,which.max)
    
    #test group assignment and y prediction
    pai_est=sapply(1:K, function(k) 
      exp(cbind(1,G_test) %*% gamma_est_matrix[,k,drop=F])/(rowSums(exp(cbind(1,G_test) %*% gamma_est_matrix))))
    Y_prd_test=apply(sapply(1:K, function(x) pai_est[,x]*(beta0_est[x]+X_test %*% beta_est)),1,sum)
    
    f_est<-sapply(1:K, function(x) (1/sqrt(2*pi*sigma2_est))*exp(-(Y_test-beta0_est[x]-X_test %*% beta_est)^2/(2*sigma2_est)))
    #calculate the expected value of Z
    w_est=sapply(1:K, function(k) 
      (pai_est[,k]*f_est[,k])/diag(pai_est %*% t(f_est)))
    grp_assign_test<-apply(w_est,1,which.max)
    #Y_prd_test <- beta0_est[grp_assign_test] + X_test %*% beta_est
    
    onecv.res<-list(K=best.K, lambda=best.lambda, fld_id=flds[[f]], grp_assign_train=grp_assign_train,
                    grp_assign_test=grp_assign_test, Y_prd_test=Y_prd_test, BIC=BIC.all, theta_est=best.par)
    cvres.our1000[[f]]<-onecv.res
  }
  cvres.our1000.all[[ii]]<-cvres.our1000
}
save.image(file="lung_disease_our_cv_top2000_fast_K3_robust_compare2.Rdata")

n=nrow(gdata2)

Y_prd<-Y
grp_assign<-Y
RMSE<-c()
R2<-c()
for(ii in 1:length(lambda.all)){
  for(i in 1:length(flds)){
    id=cvres.our1000.all[[ii]][[i]]$fld_id
    Y_prd[id]<-cvres.our1000.all[[ii]][[i]]$Y_prd_test
    grp_assign[id]<-cvres.our1000.all[[ii]][[i]]$grp_assign_test
  }
  RMSE[ii]=sqrt(sum((Y_prd-Y)^2)/n)
  R2[ii]=1-sum((Y-Y_prd)^2)/sum((Y-mean(Y))^2)
}
RMSE
R2

#load("D:/research/Test of depenence with Kayhan and mingming/Lung disease data/diagnosis_iPF_paper.Rdata")
grp_diag<-ifelse(LungData.all$diagnosis$diagnosis=="COPD",1,2)[index.Y]
table(grp_diag, grp_assign)

pvalues=kruskal.test(Y ~ grp_assign)$p.value
pvalues
