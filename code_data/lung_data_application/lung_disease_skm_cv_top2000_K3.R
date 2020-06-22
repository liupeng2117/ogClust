##lung disease##
load("/home/pel67/subtyping/lung/data/LungData.all.Rdata")
library(KLClust)
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

cvres.skm<-foreach(f=1:length(flds),.packages = "KLClust") %dopar% {
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

  
  #choose the optimal k and lambda
  kl<-KL.S4(t(gdata2_train),lambda_list = list(seq(4,15,1)), 
            k_vector=3,num.cores = 1)
  #skm clustering
  cl<-KMeansSparseCluster(gdata2_train,K=kl$optimal_k,wbounds=kl$optimal_lambda)
  weight<-cl[[1]]$ws #gene weight
  grp_assign_train<-cl[[1]]$Cs #group assignement
  n.k=kl$optimal_k
  
  #test the association between cluster assignment and clinical variables
  pvalues<-c()
  for(i in 1:ncol(LungData.all$clinic)){
    values=LungData.all$clinic[index.Y,i][-flds[[f]]]
    tempdata=cbind(clusters=grp_assign_train,values=values)
    pvalues[i]<-kruskal.test(values ~ clusters, data = tempdata)$p.value 
  }
  names(pvalues)<-colnames(LungData.all$clinic)
  
  #distance matrix
  dist.matrix=matrix(,nrow=nrow(gdata2_test),ncol=n.k)
  for(k in 1: n.k){
    miu_k<-apply(gdata2_train[grp_assign_train==k,],2,mean)
    dist.matrix[,k]<-apply(gdata2_test,1,function(x) weight %*% (x-miu_k)^2)
  }
  
  #group assignment of test data
  grp_assign_test=apply(dist.matrix, 1, which.min)
  #grp_assign_prd[flds[[f]]]=grp_assign_test
  
  #calculate the y_est
  Y_prd_test=Y_test
  for(k in 1: n.k){
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
  onecv.res<-list(pvalues=pvalues,K=n.k,lambda=kl$optimal_lambda, fld_id=flds[[f]], 
                  grp_assign=grp_assign_test, Y_prd_test=Y_prd_test)
  return(onecv.res)
}
save.image(file="lung_disease_skm_cv_top2000_K3.Rdata")

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