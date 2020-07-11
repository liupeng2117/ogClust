#---------parameters------------
setwd("D:/research/Disease subtyping/simulation/suvival_yusi/new")
n.rep=100
NG=1000
n=600 # sample
np=2 # 2 covariates
beta1=1 # coefficient of clinical var1
beta2=1 # coefficient of clinical var2
gamma1=c(3,3,3,3,3,0,0,0,0,0,-3,-3,-3,-3,-3) # 1st set coefficient of genes 1 to 15, others are 0
gamma2=c(0,0,0,0,0,3,3,3,3,3,-3,-3,-3,-3,-3) # 2nd set coefficient of genes 1 to 15, others are 0
miu=c(1,4,7) # overall mean
sigma2=1 # overall variance

# ------------------------------
# ------- ogclust --------
#-------------------------------
load(file = "survival_ogclust_K3_delta1_gamma1.rdata")

K.all=3 #Set K = truth
lambda.all<-c(0.01,0.03,0.05,0.08,0.10,0.12,0.14,0.16,0.18,0.2)
tune_pars=cbind(K=rep(K.all,each=length(lambda.all)),lambda=rep(lambda.all,length(K.all)))

#---------------------
n_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
n_true_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2 <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)

for(j in 1:nrow(tune_pars)){
  for (i in 1:n.rep) {
    #train
    Y <- gen.train$data[[i]]$time
    Z <- gen.train$data[[i]]$Z
    X1 <- gen.train$data[[i]]$X1
    X2 <- gen.train$data[[i]]$X2
    #Y_prd <- par.est.all[[i]][[j]]$time_prd
    #e=Y-(miu[Z]+beta1*X1+beta2*X2)
    #RMSE[j,i]<- sqrt(sum((Y_prd - (Y-e))^2) / n)
    #R2[j,i] <- 1 - sum((Y_prd-(Y-e))^2) / sum(((Y-e) - mean(Y))^2)
    ARI[j,i] <- par.est.all[[i]][[j]]$ARI
    n_markers[j,i]<- par.est.all[[i]][[j]]$n_genes_selected
    n_true_markers[j,i]<- par.est.all[[i]][[j]]$n_true_genes_selected 
    grp_assign_train <- par.est.all[[i]][[j]]$grp_assign_train
    # pvalues[i]=kruskal.test(Y ~ grp_assign_train)$p.value
    print(table(grp_assign_train,gen.train$grp_assign[i,]))
    
    #test
    Y_test <- gen.test$data[[i]]$time
    Z_test <- gen.test$data[[i]]$Z
    X1_test <- gen.test$data[[i]]$X1
    X2_test <- gen.test$data[[i]]$X2
    Y_prd_test <- par.est.all[[i]][[j]]$Y_prd_test
    #e_test=Y_test-(miu[Z_test]+beta1*X1_test+beta2*X2_test)
    RMSE_test[j,i]<- sqrt(sum((Y_prd_test - Y_test)^2) / n)
    R2_test[j,i] <- 1 - sum((Y_prd_test-Y_test)^2) / sum((Y_test - mean(Y_test))^2)
    ARI_test[j,i] <- par.est.all[[i]][[j]]$ARI_test
  }
}

n_markers
n_true_markers

df.ogclust<-data.frame(group=rep(1:n.rep,each=nrow(tune_pars)),
                     n_markers=as.numeric(n_markers),
                     n_true_markers=as.numeric(n_true_markers),
                     FPs=as.numeric(n_markers)-as.numeric(n_true_markers),
                     FNs=15-as.numeric(n_true_markers),
                     perc_true_markers=as.numeric(n_true_markers)/15,
                     specificity=(985-(as.numeric(n_markers)-as.numeric(n_true_markers)))/985,
                     ARI=as.numeric(ARI),
                     ARI_test=as.numeric(ARI_test),
                     R2=as.numeric(R2),
                     RMSE=as.numeric(RMSE),
                     R2_test=as.numeric(R2_test),
                     RMSE_test=as.numeric(RMSE_test)
                     )
df.ogclust$youden<-with(df.ogclust, perc_true_markers+specificity-1)
df.ogclust<-df.ogclust[df.ogclust$n_markers>10,]
library(ggplot2)
p1 <- ggplot(df.ogclust, aes(n_markers, ARI, group = group)) + geom_line() + geom_smooth(aes(group=1))
p1
p2 <- ggplot(df.ogclust, aes(n_markers, RMSE_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p2
p3 <- ggplot(df.ogclust, aes(n_markers, R2_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p3
p4 <- ggplot(df.ogclust, aes(n_markers, youden, group = group)) + geom_line() + geom_smooth(aes(group=1))
p4

x_ogclust <- df.ogclust$n_markers
y_ogclust <- df.ogclust$ARI
y2_ogclust <- df.ogclust$RMSE_test
y3_ogclust <- df.ogclust$R2_test
#y4_ogclust <- df.ogclust$youden
y4_ogclust <- df.ogclust$FPs
y5_ogclust <- df.ogclust$FNs
lo <- loess(y_ogclust~x_ogclust)
y_prd_ogclust<-predict(lo,seq(10,150,5))
lo <- loess(y2_ogclust~x_ogclust)
y2_prd_ogclust<<-predict(lo,seq(10,150,5))
lo <- loess(y3_ogclust~x_ogclust)
y3_prd_ogclust<<-predict(lo,seq(10,150,5))
lo <- loess(y4_ogclust~x_ogclust)
y4_prd_ogclust<<-predict(lo,seq(10,150,5))
lo <- loess(y5_ogclust~x_ogclust)
y5_prd_ogclust<<-predict(lo,seq(10,150,5))
# ------------------------------
# ------- skmeans --------
#-------------------------------
load(file = "survival_skmeans_K3_delta1_gamma1.rdata")

K.all=3 #Set K = truth
lambda.all<-seq(1.1,3,0.3) # penalty parameter
tune_pars=cbind(K=rep(K.all,each=length(lambda.all)),lambda=rep(lambda.all,length(K.all)))

n_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
n_true_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2 <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)

for(j in 1:nrow(tune_pars)){
  for (i in 1:n.rep) {
    #train
    Y <- gen.train$data[[i]]$time
    Z <- gen.train$data[[i]]$Z
    X1 <- gen.train$data[[i]]$X1
    X2 <- gen.train$data[[i]]$X2
    ARI[j,i] <- res.skm[[i]][[j]]$ARI
    n_markers[j,i]<- res.skm[[i]][[j]]$n_genes_selected
    n_true_markers[j,i]<- res.skm[[i]][[j]]$n_true_genes_selected 
    grp_assign_train <- res.skm[[i]][[j]]$grp_assign_train
    # pvalues[i]=kruskal.test(Y ~ grp_assign_train)$p.value
    print(table(grp_assign_train,gen.train$grp_assign[i,]))
    
    #test
    Y_test <- gen.test$data[[i]]$time
    Z_test <- gen.test$data[[i]]$Z
    X1_test <- gen.test$data[[i]]$X1
    X2_test <- gen.test$data[[i]]$X2
    Y_prd_test <- res.skm[[i]][[j]]$Y_prd_test
    #e_test=Y_test-(miu[Z_test]+beta1*X1_test+beta2*X2_test)
    RMSE_test[j,i]<- sqrt(sum((Y_prd_test - Y_test)^2) / n)
    R2_test[j,i] <- 1 - sum((Y_prd_test-Y_test)^2) / sum((Y_test - mean(Y_test))^2)
    ARI_test[j,i] <- res.skm[[i]][[j]]$ARI_test
  }
}

df.skm<-data.frame(group=rep(1:n.rep,each=nrow(tune_pars)),
                     n_markers=as.numeric(n_markers),
                     n_true_markers=as.numeric(n_true_markers),
                     FPs=as.numeric(n_markers)-as.numeric(n_true_markers),
                     FNs=15-as.numeric(n_true_markers),
                     perc_true_markers=as.numeric(n_true_markers)/15,
                     specificity=(985-(as.numeric(n_markers)-as.numeric(n_true_markers)))/985,
                     ARI=as.numeric(ARI),
                     ARI_test=as.numeric(ARI_test),
                     R2=as.numeric(R2),
                     RMSE=as.numeric(RMSE),
                     R2_test=as.numeric(R2_test),
                     RMSE_test=as.numeric(RMSE_test)
)
df.skm$youden<-with(df.skm, perc_true_markers+specificity-1)

library(ggplot2)
p1 <- ggplot(df.skm, aes(n_markers, ARI, group = group)) + geom_line() + geom_smooth(aes(group=1))
p1
p2 <- ggplot(df.skm, aes(n_markers, RMSE_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p2
p3 <- ggplot(df.skm, aes(n_markers, R2_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p3
p4 <- ggplot(df.skm, aes(n_markers, youden, group = group)) + geom_line() + geom_smooth(aes(group=1))
p4

x_skm <- df.skm$n_markers
y_skm <- df.skm$ARI
y2_skm <- df.skm$RMSE_test
y3_skm <- df.skm$R2_test
#y4_skm <- df.skm$youden
y4_skm <- df.skm$FPs
y5_skm <- df.skm$FNs
lo <- loess(y_skm~x_skm)
y_prd_skm<-predict(lo,seq(10,150,5))
lo <- loess(y2_skm~x_skm)
y2_prd_skm<<-predict(lo,seq(10,150,5))
lo <- loess(y3_skm~x_skm)
y3_prd_skm<<-predict(lo,seq(10,150,5))
lo <- loess(y4_skm~x_skm)
y4_prd_skm<<-predict(lo,seq(10,150,5))
lo <- loess(y5_skm~x_skm)
y5_prd_skm<<-predict(lo,seq(10,150,5))

# ------------------------------
# ------- pmbc --------
#-------------------------------
load(file = "survival_pmbc_K3_delta1_gamma1.rdata")

K.all=3 #Set K = truth
lambda.all<-seq(16,24,2) # penalty parameter
tune_pars=cbind(K=rep(K.all,each=length(lambda.all)),lambda=rep(lambda.all,length(K.all)))

n_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
n_true_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2 <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)

for(j in 1:nrow(tune_pars)){
  for (i in 1:n.rep) {
    #train
    Y <- gen.train$data[[i]]$time
    Z <- gen.train$data[[i]]$Z
    X1 <- gen.train$data[[i]]$X1
    X2 <- gen.train$data[[i]]$X2
    ARI[j,i] <- res.pmbc[[i]][[j]]$ARI
    n_markers[j,i]<- res.pmbc[[i]][[j]]$n_genes_selected
    n_true_markers[j,i]<- res.pmbc[[i]][[j]]$n_true_genes_selected 
    grp_assign_train <- res.pmbc[[i]][[j]]$grp_assign_train
    # pvalues[i]=kruskal.test(Y ~ grp_assign_train)$p.value
    #print(table(grp_assign_train,gen.train$grp_assign[i,]))
    
    #test
    Y_test <- gen.test$data[[i]]$time
    Z_test <- gen.test$data[[i]]$Z
    X1_test <- gen.test$data[[i]]$X1
    X2_test <- gen.test$data[[i]]$X2
    Y_prd_test <- res.pmbc[[i]][[j]]$Y_prd_test
    #e_test=Y_test-(miu[Z_test]+beta1*X1_test+beta2*X2_test)
    RMSE_test[j,i]<- sqrt(sum((Y_prd_test - Y_test)^2) / n)
    R2_test[j,i] <- 1 - sum((Y_prd_test-Y_test)^2) / sum((Y_test - mean(Y_test))^2)
    ARI_test[j,i] <- res.pmbc[[i]][[j]]$ARI_test
  }
}

df.pmbc<-data.frame(group=rep(1:n.rep,each=nrow(tune_pars)),
                      n_markers=as.numeric(n_markers),
                      n_true_markers=as.numeric(n_true_markers), 
                      FPs=as.numeric(n_markers)-as.numeric(n_true_markers),
                      FNs=15-as.numeric(n_true_markers),
                      perc_true_markers=as.numeric(n_true_markers)/15,
                      specificity=(985-(as.numeric(n_markers)-as.numeric(n_true_markers)))/985,
                      ARI=as.numeric(ARI),
                      ARI_test=as.numeric(ARI_test),
                      R2=as.numeric(R2),
                      RMSE=as.numeric(RMSE),
                      R2_test=as.numeric(R2_test),
                      RMSE_test=as.numeric(RMSE_test)
)
df.pmbc$youden<-with(df.pmbc, perc_true_markers+specificity-1)

library(ggplot2)
p1 <- ggplot(df.pmbc, aes(n_markers, ARI, group = group)) + geom_line() + geom_smooth(aes(group=1))
p1
p2 <- ggplot(df.pmbc, aes(n_markers, RMSE_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p2
p3 <- ggplot(df.pmbc, aes(n_markers, R2_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p3
p4 <- ggplot(df.pmbc, aes(n_markers, youden, group = group)) + geom_line() + geom_smooth(aes(group=1))
p4

x_pmbc <- df.pmbc$n_markers
y_pmbc <- df.pmbc$ARI
y2_pmbc <- df.pmbc$RMSE_test
y3_pmbc <- df.pmbc$R2_test
#y4_pmbc <- df.pmbc$youden
y4_pmbc <- df.pmbc$FPs
y5_pmbc <- df.pmbc$FNs
lo <- loess(y_pmbc~x_pmbc)
y_prd_pmbc<-predict(lo,seq(10,150,5))
lo <- loess(y2_pmbc~x_pmbc)
y2_prd_pmbc<<-predict(lo,seq(10,150,5))
lo <- loess(y3_pmbc~x_pmbc)
y3_prd_pmbc<<-predict(lo,seq(10,150,5))
lo <- loess(y4_pmbc~x_pmbc)
y4_prd_pmbc<<-predict(lo,seq(10,150,5))
lo <- loess(y5_pmbc~x_pmbc)
y5_prd_pmbc<<-predict(lo,seq(10,150,5))

# ------------------------------
# ------- scluster --------
#-------------------------------
load(file = "survival_scluster_K3_delta1_gamma1.rdata")

K.all=3 #Set K = truth
ngs<-seq(5,200,10) # penalty parameter
tune_pars=cbind(K=rep(K.all,each=length(ngs)),ngs=rep(ngs,length(K.all)))

n_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
n_true_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2 <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)

for(j in 1:nrow(tune_pars)){
  for (i in 1:n.rep) {
    if(i!=4){
      #train
      Y <- gen.train$data[[i]]$time
      Z <- gen.train$data[[i]]$Z
      X1 <- gen.train$data[[i]]$X1
      X2 <- gen.train$data[[i]]$X2
      ARI[j,i] <- res.scluster[[i]][[j]]$ARI
      n_markers[j,i]<- res.scluster[[i]][[j]]$n_genes_selected
      n_true_markers[j,i]<- res.scluster[[i]][[j]]$n_true_genes_selected 
      grp_assign_train <- res.scluster[[i]][[j]]$grp_assign_prd
      # pvalues[i]=kruskal.test(Y ~ grp_assign_train)$p.value
      #print(table(grp_assign_train,gen.train$grp_assign[i,]))
      
      #test
      Y_test <- gen.test$data[[i]]$time
      Z_test <- gen.test$data[[i]]$Z
      X1_test <- gen.test$data[[i]]$X1
      X2_test <- gen.test$data[[i]]$X2
      Y_prd_test <- res.scluster[[i]][[j]]$Y_prd_test
      #e_test=Y_test-(miu[Z_test]+beta1*X1_test+beta2*X2_test)
      RMSE_test[j,i]<- sqrt(sum((Y_prd_test - Y_test)^2) / n)
      R2_test[j,i] <- 1 - sum((Y_prd_test-Y_test)^2) / sum((Y_test - mean(Y_test))^2)
      ARI_test[j,i] <- res.scluster[[i]][[j]]$ARI_test
    }
  }
}

df.scluster<-data.frame(group=rep(1:(n.rep-1),each=nrow(tune_pars)),
                        n_markers=as.numeric(n_markers[,-4]),
                        n_true_markers=as.numeric(n_true_markers[,-4]), 
                        FPs=as.numeric(n_markers[,-4])-as.numeric(n_true_markers[,-4]),
                        FNs=15-as.numeric(n_true_markers[,-4]),
                        perc_true_markers=as.numeric(n_true_markers[,-4])/15,
                        specificity=(985-(as.numeric(n_markers[,-4])-as.numeric(n_true_markers[,-4])))/985,
                        ARI=as.numeric(ARI[,-4]),
                        ARI_test=as.numeric(ARI_test[,-4]),
                        R2=as.numeric(R2[,-4]),
                        RMSE=as.numeric(RMSE[,-4]),
                        R2_test=as.numeric(R2_test[,-4]),
                        RMSE_test=as.numeric(RMSE_test[,-4])
                       )

df.scluster$youden<-with(df.scluster, perc_true_markers+specificity-1)

library(ggplot2)
p1 <- ggplot(df.scluster, aes(n_markers, ARI_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p1
p2 <- ggplot(df.scluster, aes(n_markers, RMSE_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p2
p3 <- ggplot(df.scluster, aes(n_markers, R2_test, group = group)) + geom_line() + geom_smooth(aes(group=1))
p3
p4 <- ggplot(df.scluster, aes(n_markers, youden, group = group)) + geom_line() + geom_smooth(aes(group=1))
p4

x_scluster <- df.scluster$n_markers
y_scluster <- df.scluster$ARI
y2_scluster <- df.scluster$RMSE_test
y3_scluster <- df.scluster$R2_test
#y4_scluster <- df.scluster$youden
y4_scluster <- df.scluster$FPs
y5_scluster <- df.scluster$FNs
lo <- loess(y_scluster~x_scluster)
y_prd_scluster<-predict(lo,seq(10,150,5))
lo <- loess(y2_scluster~x_scluster)
y2_prd_scluster<<-predict(lo,seq(10,150,5))
lo <- loess(y3_scluster~x_scluster)
y3_prd_scluster<<-predict(lo,seq(10,150,5))
lo <- loess(y4_scluster~x_scluster)
y4_prd_scluster<<-predict(lo,seq(10,150,5))
lo <- loess(y5_scluster~x_scluster)
y5_prd_scluster<<-predict(lo,seq(10,150,5))

#id_scluster<-order(x_scluster)
#id_ogclust<-order(x_ogclust)
#id_skm<-order(x_skm)
#id_pmbc<-order(x_pmbc)

df1=rbind(data.frame(x=seq(10,150,5),y=y_prd_ogclust,group="ogClust"),
          data.frame(x=seq(10,150,5),y=y_prd_scluster,group="SC"),
          data.frame(x=seq(10,150,5),y=y_prd_skm,group="SKM"),
          data.frame(x=seq(10,150,5),y=y_prd_pmbc,group="PMBC"))
p1<-ggplot(df1,aes(x,y,group)) +
  geom_line(aes(color=group, linetype=group),size=1)+
  #geom_point(aes(color=group, shape=group))+
  xlab("number of markers") + ylab("ARI") +
  ylim(0,1)+
  theme_classic()

df2=rbind(data.frame(x=seq(10,150,5),y=y2_prd_ogclust,group="ogClust"),
          data.frame(x=seq(10,150,5),y=y2_prd_scluster,group="SC"),
          data.frame(x=seq(10,150,5),y=y2_prd_skm,group="SKM"),
          data.frame(x=seq(10,150,5),y=y2_prd_pmbc,group="PMBC"))
p2<-ggplot(df2,aes(x,y,group)) +
  geom_line(aes(color=group, linetype=group), size=1)+
  #geom_point(aes(color=group, shape=group))+
  xlab("number of markers") + ylab("RMSE") +
  ylim(0,30)+
  theme_classic()

df3=rbind(data.frame(x=seq(10,150,5),y=y3_prd_ogclust,group="ogClust"),
          data.frame(x=seq(10,150,5),y=y3_prd_scluster,group="SC"),
          data.frame(x=seq(10,150,5),y=y3_prd_skm,group="SKM"),
          data.frame(x=seq(10,150,5),y=y3_prd_pmbc,group="PMBC"))
p3<-ggplot(df3,aes(x,y,group)) +
  geom_line(aes(color=group,linetype=group), size=1)+
  #geom_point(aes(color=group, shape=group))+
  xlab("number of markers") + ylab("R2") +
  ylim(0,1)+
  theme_classic()

df4=rbind(data.frame(x=seq(10,150,5),y=y4_prd_ogclust,group="ogClust"),
          data.frame(x=seq(10,150,5),y=y4_prd_scluster,group="SC"),
          data.frame(x=seq(10,150,5),y=y4_prd_skm,group="SKM"),
          data.frame(x=seq(10,150,5),y=y4_prd_pmbc,group="PMBC"))
p4<-ggplot(df4,aes(x,y,group)) +
  geom_line(aes(color=group,linetype=group), size=1)+
  #geom_point(aes(color=group, shape=group))+
  xlab("number of markers") + ylab("FPs") +
  ylim(0,150)+
  theme_classic()

df5=rbind(data.frame(x=seq(10,150,5),y=y5_prd_ogclust,group="ogClust"),
          data.frame(x=seq(10,150,5),y=y5_prd_scluster,group="SC"),
          data.frame(x=seq(10,150,5),y=y5_prd_skm,group="SKM"),
          data.frame(x=seq(10,150,5),y=y5_prd_pmbc,group="PMBC"))
p5<-ggplot(df5,aes(x,y,group)) +
  geom_line(aes(color=group,linetype=group), size=1)+
  #geom_point(aes(color=group, shape=group))+
  xlab("number of markers") + ylab("FNs") +
  ylim(0,15)+
  theme_classic()

df.delta1.gamma1<-rbind(data.frame(measures="ARI",df1),
                 data.frame(measures="RMSE",df2),
                 data.frame(measures="R2",df3),
                 data.frame(measures="FPs",df4),
                 data.frame(measures="FNs",df5))
save(df.delta1.gamma1, file="df.delta1.gamma1.rdata")
#library(gridExtra)
#install.packages("ggpubr")
library(ggpubr)
#grid.arrange(p1,p2,p3,p4,nrow=2)
pdf(file="smooth.plot.normal2.pdf",width = 12, height=9)
ggarrange(p1,p2,p3,p4,p5,labels=c("A","B","C","D","E"),nrow=3,ncol=2)
dev.off()
