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
gdata2=gdata[index.Y,]
dim(gdata2) #234 x 2000

#covariates to adjust for 
X1=age=LungData.all$clinic[index.Y,"AGE"] 
X2=bmi=LungData.all$clinic[index.Y,"BMI"] 
X3=gender=LungData.all$clinic[index.Y,"GENDER"] #gender
X=cbind(X1,X2,X3)
# number of genes and patients
n=nrow(gdata2)
NG=ncol(gdata2)
np=ncol(X)

library(doParallel) # foreach
cl<-makeCluster(40)
registerDoParallel(cl)

B=500
boot.res<-foreach(b=1:B,.packages="KLClust") %dopar% {
  print(paste0("b=",b))
  id<-sample(1:n,n,replace = T)
  age.b=age[id]
  bmi.b=bmi[id]
  gender.b=gender[id]
  gdata2.b=gdata2[id,]

  #choose the optimal k and lambda
  kl<-KL.S4(t(gdata2.b),lambda_list = list(seq(4,10,3)), 
            k_vector=c(2),num.cores = 1)
  #skm clustering
  cl<-KMeansSparseCluster(gdata2.b,K=kl$optimal_k,wbounds=kl$optimal_lambda)
  weight<-cl[[1]]$ws #gene weight
  grp_assign<-cl[[1]]$Cs #group assignement
  n.k=kl$optimal_k

  final.res<-list(K=n.k, grp_assign=grp_assign, weight=weight)
  return(final.res)
}
save.image(file="lung_disease_skm_top2000_bootstrap_B500.Rdata")