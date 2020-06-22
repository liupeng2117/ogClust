##lung disease##
source("/home/pel67/subtyping/lung/code/EM_function_fast_elasticnet.R")
source("/home/pel67/subtyping/lung/code/fit.mixture.model_function_elasticnet.R")
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
K.all=3
#lambda.all=c(2,4,6,8,10)
lambda.all=c(0.187637400, 0.13, 0.11 ,0.097834296, 0.07, 0.051010883, 0.026597117)
tune_pars=cbind(K=rep(K.all,each=length(lambda.all)),lambda=rep(lambda.all,length(K.all)))
# number of genes and patients
n=nrow(G)
NG=ncol(G)
np=ncol(X)
  
res.all<-foreach(l=1:nrow(tune_pars)) %dopar% {
    #parameters
    K=tune_pars[l,"K"]
    lambda=tune_pars[l,"lambda"]
    #print(paste0("K=",K))
    #res.oneK<-foreach(r=1:length(lambda.all),.combine=rbind) %dopar% {
    #lambda=lambda.all[r]
    #2 initials and choose the best
    n.int=2
    est.par=matrix(,ncol=((K-1)*(NG+1)+np+K+6),nrow=n.int)
    #colnames(est.par)<-colnames(par.all)
    theta.int.all=matrix(,ncol=(K-1)*(NG+1)+np+K+1,nrow=n.int)
    for(i in 1:n.int){
        fit.res<-fit.mixture.model(n=n, K=K, np=np, NG=NG, lambda=lambda, 
                                   alpha=0.5, df=df, G=G, Y=Y, X=X)
        est.par[i,]<-fit.res$res
    }
    #choose the best par estimates which has the maximum ll
    best.par<-est.par[which.min(est.par[,(K-1)*(NG+1)+K+np+5]),]
    #best.id<-which.min(best.par[,(K-1)*(NG+1)+K+np+5])
    #res.all[[l]]<-best.par[best.id,]
    return(best.par)
}

print("All K and lambda DONE")
BIC.all<-c()
for (l in 1:nrow(tune_pars)){
    K=tune_pars[l,"K"]
    BIC.all[l]<-res.all[[l]][(K-1)*(NG+1)+np+K+5]
} 

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
  
#group assignment and y prediction
pai_est=sapply(1:K, function(k) 
    exp(cbind(1,G) %*% gamma_est_matrix[,k,drop=F])/(rowSums(exp(cbind(1,G) %*% gamma_est_matrix))))
f_est<-sapply(1:K, function(x) (1/sqrt(2*pi*sigma2_est))*exp(-(Y-beta0_est[x]-X %*% beta_est)^2/(2*sigma2_est)))
Y_prd=apply(sapply(1:K, function(x) pai_est[,x]*(beta0_est[x]+X %*% beta_est)),1,sum)
#calculate the expected value of Z
w_est=sapply(1:K, function(k) 
    (pai_est[,k]*f_est[,k])/diag(pai_est %*% t(f_est)))
grp_assign<-apply(pai_est,1,which.max)

#pvalues
tempdata=cbind(clusters=grp_assign,values=Y)
pvalues=kruskal.test(values ~ clusters, data = tempdata)$p.value
  
final.res<-list(pvalues=pvalues, K=best.K, lambda=best.lambda, grp_assign=grp_assign,
                Y_prd=Y_prd, BIC=BIC.all, theta_est=best.par)
#-------
save.image(file="lung_disease_our_top2000_elastic_fast_K3.Rdata")

sqrt(sum((Y_prd - Y)^2)/n)
1-sum((Y - Y_prd)^2)/sum((Y - mean(Y))^2)

grp_diag<-ifelse(LungData.all$diagnosis$diagnosis=="COPD",1,2)[index.Y]
table(grp_diag, grp_assign)
     
kruskal.test(Y ~ grp_assign)$p.value

fisher.pvalue=fisher.test(factor(grp_diag), factor(grp_assign))$p.value
fisher.pvalue

#load lung_disease_our_top2000_elastic_fast_K3.Rdata, set best.par=res.all[[5]]
plot(pai_est[,1],pai_est[,3],col=grp_diag, xlab="p1", ylab="p3")
legend("topright",legend=c("COPD","ILD"), col=c(1,2),pch = 1)
color<-grp_assign
color[grp_assign==2]<-1
color[grp_assign==1]<-2
plot(pai_est[,1],pai_est[,3],col=color,xlab="p1",ylab="p3")
legend("topright",legend=c("group1","group2","group3"), col=c(2,1,3),pch = 1)

df[,-1]<-apply(df[,-1],2,function(x) (x-mean(x))/sd(x))
#pca
pca <- prcomp(pai_est,
              center = TRUE,
              scale. = TRUE) 
dataset=data.frame(pca$x)
prop.pca = pca$sdev^2/sum(pca$sdev^2)

#predict
group=as.factor(grp_assign)
require(scales)
ggplot(dataset) + geom_point(aes(PC1, PC2, colour = group, shape = group), size = 2.5) +
  scale_color_manual(values=c("red", "black", "green"))+
  labs(x = paste("PC1 (", percent(prop.pca[1]), ")", sep=""),
       y = paste("PC2 (", percent(prop.pca[2]), ")", sep="")) +
  theme_bw()

#diagnosis
group=as.factor(grp_diag)
require(scales)
ggplot(dataset) + geom_point(aes(PC1, PC2, colour = group, shape = group), size = 2.5) +
  scale_color_manual(values=c("black","red"))+
  labs(x = paste("PC1 (", percent(prop.pca[1]), ")", sep=""),
       y = paste("PC2 (", percent(prop.pca[2]), ")", sep="")) +
  theme_bw()

# dot plot
grp1.index<-which(grp_assign==1)
grp2.index<-which(grp_assign==2)
grp3.index<-which(grp_assign==3)
data<-cbind(pai_est,grp_assign)
data2<-data[c(grp1.index,grp2.index,grp3.index),]
data2<-cbind(1:length(grp_assign),data2)
colnames(data2)<-c("index","p1","p2","p3","grp")
data2<-data.frame(data2)
library(reshape2)
data2.long<-melt(data2,id.vars=c("grp","index"),variable.names="prob_grp",value.name="probability")

color<-numeric(length(data2.long$variable))
color[data2.long$variable=="p1"]<-2
color[data2.long$variable=="p2"]<-1
color[data2.long$variable=="p3"]<-3
with(data2.long, plot(index,probability,col=color))

#boxplot
grp_assign2<-mapvalues(grp_assign,from=c(1,2,3),to=c(2,1,3))
df=data.frame(grp_assign2,Y)
ggplot(df,aes(x=grp_assign2,y=Y,group=grp_assign2)) + geom_boxplot()+xlab("Group")+ylab("FEV1%prd")+theme_classic()

