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

n=nrow(gdata2)
#choose the optimal k and lambda
kl<-KL.S4(t(gdata2),lambda_list = list(seq(4,15,1)), 
            k_vector=3,num.cores = 1)
#skm clustering
cl<-KMeansSparseCluster(gdata2,K=kl$optimal_k,wbounds=kl$optimal_lambda)
weight<-cl[[1]]$ws #gene weight
grp_assign<-cl[[1]]$Cs #group assignement
n.k=kl$optimal_k
  
#test the association between cluster assignment and clinical variables
pvalues<-c()
for(i in 1:ncol(LungData.all$clinic)){
    values=LungData.all$clinic[index.Y,i]
    tempdata=cbind(clusters=grp_assign,values=values)
    pvalues[i]<-kruskal.test(values ~ clusters, data = tempdata)$p.value 
}
names(pvalues)<-colnames(LungData.all$clinic)

res<-list(pvalues=pvalues,K=n.k,lambda=kl$optimal_lambda,
                  grp_assign=grp_assign, Y_prd=Y_prd)

save.image(file="lung_disease_skm_top2000_K3.Rdata")

weight

grp_diag<-ifelse(LungData.all$diagnosis$diagnosis=="COPD",1,2)[index.Y]
table(grp_diag, grp_assign)

pvalues=kruskal.test(Y ~ grp_assign)$p.value
pvalues
#1.083167e-06
fisher.pvalue=fisher.test(factor(grp_diag), factor(grp_assign))$p.value
fisher.pvalue

pvalues.all<-c()
for(i in 1:ncol(LungData.all$clinic)){
  yy<-LungData.all$clinic[index.Y,i]
  pvalues.all[i]=kruskal.test(yy ~ grp_assign)$p.value
}

names(pvalues.all)<-colnames(LungData.all$clinic)
pvalues.all["GENDER"]<-chisq.test(as.factor(LungData.all$clinic$GENDER[index.Y]), as.factor(grp_assign))$p.value
pvalues.all[c("AGE_yrs","pkyrs","BMI","GENDER")]
sort(pvalues.all)[1:10]

pvalues.all2<-c()
for(i in c(1:6)){
  yy <- LungData.all$diagnosis[index.Y,i]
  pvalues.all2[i]=chisq.test(as.factor(yy), as.factor(grp_assign))$p.value
}

#boxplot
grp_assign2<-mapvalues(grp_assign,from=c(1,2,3),to=c(2,3,1))
df=data.frame(grp_assign2=factor(grp_assign2),Y)
ggplot(df,aes(x=grp_assign2,y=Y,group=grp_assign2, fill=grp_assign2)) + geom_boxplot()+xlab("Group")+ylab("FEV1%prd")+theme_classic()


