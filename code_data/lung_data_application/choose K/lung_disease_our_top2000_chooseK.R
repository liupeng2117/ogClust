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

# Y: exlude patients with missing outcome
index.Y=!is.na(LungData.all$clinic[,"fev1pd1a"]) #exclude if outcome is missing
Y=LungData.all$clinic[index.Y,"fev1pd1a"]
gdata2=gdata[index.Y,]
dim(gdata2) #234 x 2000
G=gdata2

#covariates X to adjust for 
X1=age=LungData.all$clinic[index.Y,"AGE"] 
X2=bmi=LungData.all$clinic[index.Y,"BMI"] 
X3=gender=LungData.all$clinic[index.Y,"GENDER"] #gender
X=cbind(age,bmi, gender)
#data: X Y G
df=data.frame(X1=age,X2=bmi,X3=gender,G=G,Y=Y)

#-----------------
#parallel computing
library(doParallel) # foreach
cl<-makeCluster(21)
registerDoParallel(cl)

#tuning parameters
K.all=c(1,2,3,4,5,6,7)
#lambda.all=c(2,4,6,8,10)
lambda.all=c(0.187637400, 0.13, 0.11 ,0.097834296, 0.07, 0.051010883, 0.026597117, 0.01)
tune_pars=cbind(K=rep(K.all,each=length(lambda.all)),lambda=rep(lambda.all,length(K.all)))
# number of genes and patients
n=nrow(G)
NG=ncol(G)
np=ncol(X)
  
res.all<-foreach(l=1:nrow(tune_pars),.packages = "OGClust") %dopar% {
    #parameters
    K=tune_pars[l,"K"]
    lambda=tune_pars[l,"lambda"]
    if(K==1){
      best.par=NULL
      fit=lm(Y ~ X1 + X2 + X3)
      Y_prd <- cbind(1,X) %*% fit$coefficients
      RMSE = sqrt(sum(( Y_prd - Y)^2)/n)
      R2 = 1 - sum((Y-Y_prd)^2)/sum((Y - mean(Y))^2)
      return(list(best.par=best.par, RMSE=RMSE, R2=R2, RMSE.hard=RMSE, R2.hard=R2))
    } else{
      #print(paste0("K=",K))
      #res.oneK<-foreach(r=1:length(lambda.all),.combine=rbind) %dopar% {
      #lambda=lambda.all[r]
      #2 initials and choose the best
      n.int=2
      est.par=matrix(,ncol=((K-1)*(NG+1)+np+K+6),nrow=n.int)
      #colnames(est.par)<-colnames(par.all)
      #theta.int.all=matrix(,ncol=(K-1)*(NG+1)+np+K+1,nrow=n.int)
      for(i in 1:n.int){
        beta_int <- runif(np, 0, 3)
        gamma_int <- runif((K - 1) * (NG + 1), 0, 1)
        beta0_int <- runif(K, 0, 3)
        sigma2_int <- runif(1, 1, 3)
        theta_int <- c(beta_int, gamma_int, beta0_int, sigma2_int)
        fit.res<-fit.OGClust(n=n, K=K, np=np, NG=NG, lambda=lambda, 
                             alpha=0.5, G=G, Y=Y, X=X, theta_int = theta_int,
                             robust="none")
        est.par[i,]<-fit.res$res
      }
      #choose the best par estimates which has the maximum ll
      best.par<-est.par[which.min(est.par[,(K-1)*(NG+1)+K+np+5]),]
      #choose the best par estimates 
      beta_est=best.par[1:(np)]
      gamma_est=best.par[(np+1):((K-1)*(NG+1)+np)]
      beta0_est=best.par[((K-1)*(NG+1)+np+1):((K-1)*(NG+1)+K+np)]
      sigma2_est=best.par[(K-1)*(NG+1)+np+K+1]
      gamma_est_matrix=matrix(gamma_est,ncol=K-1,byrow=TRUE)
      gamma_est_matrix=cbind(gamma_est_matrix,0)
      
      #group assignment and y prediction
      pai_est=sapply(1:K, function(k) 
        exp(cbind(1,G) %*% gamma_est_matrix[,k,drop=F])/(rowSums(exp(cbind(1,G) %*% gamma_est_matrix))))
      f_est<-sapply(1:K, function(x) (1/sqrt(2*pi*sigma2_est))*exp(-(Y-beta0_est[x]-X %*% beta_est)^2/(2*sigma2_est)))
      grp_assign<-apply(pai_est,1,which.max)
      Y_prd= beta0_est[grp_assign]+X %*% beta_est
      
      probs=sapply(1:n,function(x) pai_est[x,grp_assign[x]])
      RMSE=sqrt(sum(probs*(Y_prd - Y)^2)/sum(probs))
      R2=1-sum(probs*(Y - Y_prd)^2)/sum(probs*(Y - mean(Y))^2)
      RMSE.hard=sqrt(sum((Y_prd - Y)^2)/n)
      R2.hard=1-sum((Y - Y_prd)^2)/sum((Y - mean(Y))^2)
      return(list(best.par=best.par, RMSE=RMSE, R2=R2, RMSE.hard=RMSE.hard, R2.hard=R2.hard))
    }
}

#-------
save.image(file="lung_disease_our_top2000_chooseK.Rdata")

RMSE.all<-matrix(nrow=length(K.all),ncol=length(lambda.all))
R2.all<-matrix(nrow=length(K.all),ncol=length(lambda.all))
rownames(RMSE.all)=rownames(R2.all)=K.all
colnames(RMSE.all)=colnames(R2.all)=lambda.all
l=1
for(i in 1:length(K.all)){
  for(j in 1:length(lambda.all)){
    RMSE.all[i,j]<-res.all[[l]][["RMSE"]]
    R2.all[i,j]<-res.all[[l]][["R2"]]
    l=l+1
  }
}
RMSE.all
R2.all

RMSE.hard.all<-matrix(nrow=length(K.all),ncol=length(lambda.all))
R2.hard.all<-matrix(nrow=length(K.all),ncol=length(lambda.all))
rownames(RMSE.hard.all)=rownames(R2.all)=K.all
colnames(RMSE.hard.all)=colnames(R2.all)=lambda.all
l=1
for(i in 1:length(K.all)){
  for(j in 1:length(lambda.all)){
    RMSE.hard.all[i,j]<-res.all[[l]][["RMSE.hard"]]
    R2.hard.all[i,j]<-res.all[[l]][["R2.hard"]]
    l=l+1
  }
}
RMSE.hard.all
R2.hard.all

