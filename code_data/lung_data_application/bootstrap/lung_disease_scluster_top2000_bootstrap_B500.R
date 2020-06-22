##lung disease##
load("/home/pel67/subtyping/lung/data/LungData.all.Rdata")
library(NbClust)
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
genes<-colnames(gdata)

#standardize gene expression
gdata=sapply(1:ncol(gdata), function(x) (gdata[,x]-mean(gdata[,x]))/sd(gdata[,x]))
colnames(gdata)=genes

# exlude patients with missing outcome
# exlude genes if correlation with outcome is close to 0(<=0.2)
index.Y=!is.na(LungData.all$clinic[,"fev1pd1a"]) #exclude if outcome is missing
Y=LungData.all$clinic[index.Y,"fev1pd1a"] 
gdata2=gdata[index.Y,]
dim(gdata2) #234 x 2000

#covariates to adjust for 
X1=age=LungData.all$clinic[index.Y,"AGE"] 
X2=bmi=LungData.all$clinic[index.Y,"BMI"] 
X3=gender=LungData.all$clinic[index.Y,"GENDER"]
X=cbind(X1,X2,X3)
# number of genes and patients
n=nrow(gdata2)
NG=ncol(gdata2)
np=ncol(X)

# 5 fold cross validation
#parallel computing
library(doParallel) # foreach
cl<-makeCluster(50)
registerDoParallel(cl)

# 5 fold
library(caret)
set.seed(123)
flds <- createFolds(1:n, k = 5, list = TRUE, returnTrain = FALSE)

B=500
set.seed(123)
boot.res<-foreach(b = 1:B, .packages="NbClust") %dopar% {
  print(paste0("b=",b))
  id<-sample(1:n,n,replace = T)
  age.b=age[id]
  bmi.b=bmi[id]
  gender.b=gender[id]
  gdata2.b=gdata2[id,]
  Y.b=Y[id]

  # select the outcome related genes
  ps<-c()
  for(i in 1:ncol(gdata2.b)){
    ps[i]<-t.test(Y.b,gdata2.b[,i])$p.value
  }
  ngs<-c(5,10,20,25,30,35,40,45,50,55,60,70,80)
  RMSE<-c()
  for(o in 1:length(ngs)){
    ngsm<-ngs[o]
    gdata3.b<-gdata2.b[,order(ps)[1:ngsm]]
   cvres.scluster<-list()
   for(f in 1:length(flds)) {
      #train and test data
      gdata3_train<-gdata3.b[-flds[[f]],] #train g
      gdata3_test<-gdata3.b[flds[[f]],] #test g
      age_train<-age.b[-flds[[f]]] #train x
      age_test<-age.b[flds[[f]]] #train x
      bmi_train<-bmi.b[-flds[[f]]] #train x
      bmi_test<-bmi.b[flds[[f]]] #train x
      gender_train<-gender.b[-flds[[f]]] #train x
      gender_test<-gender.b[flds[[f]]] #train x
      Y_train=Y.b[-flds[[f]]] #train y
      Y_test=Y.b[flds[[f]]] #test y
    
    
      #kmeans clustering
      n.k<-NbClust(data=gdata3_train,min.nc = 2, max.nc = 5,method="kmeans",index="ch")$Best.nc["Number_clusters"]
      cl<-kmeans(gdata3_train,n.k)
      grp_assign_train<-cl$cluster #group assignement
    
      #distance matrix
      dist.matrix=matrix(,nrow=nrow(gdata3_test),ncol=n.k)
      for(k in 1: n.k){
        miu_k<-cl$centers[k,]
        dist.matrix[,k]<-apply(gdata3_test,1,function(x) sum((x-miu_k)^2))
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
      onecv.res<-list(K=n.k, fld_id=flds[[f]], 
                      grp_assign=grp_assign_test, Y_prd_test=Y_prd_test)
      cvres.scluster[[f]]<-onecv.res
    }
    #save.image(file="lung_disease_scluster_cv.Rdata")
  
    Y_prd<-Y.b
    grp_assign<-Y.b
    for(i in 1:length(flds)){
      id=cvres.scluster[[i]]$fld_id
      Y_prd[id]<-cvres.scluster[[i]]$Y_prd_test
      grp_assign[id]<-cvres.scluster[[i]]$grp_assign
    }
  
    RMSE[o]=sqrt(sum((Y_prd-Y)^2)/n)
    #RMSE
    #R2=1-sum((Y-Y_prd)^2)/sum((Y-mean(Y))^2)
    #R2 
  }

  final.ngs<-ngs[which.min(RMSE)]
  gdata3.b<-gdata2.b[,order(ps)[1:final.ngs]]
  n=nrow(gdata3.b)

  
  
  #kmeans clustering
  n.k<-NbClust(data=gdata3.b, min.nc = 2, max.nc = 5,method="kmeans",index="ch")$Best.nc["Number_clusters"]
  cl<-kmeans(gdata3.b,n.k)
  grp_assign<-cl$cluster #group assignement
  
  final.res<-list(K=n.k, ngs=final.ngs, grp_assign=grp_assign,genes_selected=colnames(gdata3.b))

  return(final.res)
}


save.image(file="lung_disease_scluster_top2000_bootstrap_B500.Rdata")

