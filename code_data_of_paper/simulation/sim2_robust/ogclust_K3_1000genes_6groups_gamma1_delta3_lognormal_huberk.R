rm(list=ls())
source("/home/pel67/subtyping/simulation/code/robust/new/gen.data.K3.R")
#source("~/subtyping/lung/code/EM_function_fast_elasticnet.R")
#source("~/subtyping/lung/code/fit.mixture.model_function_elasticnet.R")
library(OGClust)
# Simulation ------------------------------------------------
## Data generation --------------------
#============================
# 600 subjects, 2 subgroups, 6 clusters, 2 clinical variables x1, x2
# 1000 genes G1 to G1000
# G1 and G10 form two clusters and related to outcome
# G11 and G25 forms another three clusters but not related
# G26 to G1000 are noise
# 10% outliers on Y
#============================
set.seed(123)
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

gen.train<-gen.data.K3(n.rep, NG, n, theta, cluster_index, error="lognormal")
data.train<-gen.train$data

gen.test<-gen.data.K3(n.rep, NG, n, theta, cluster_index)
data.test<-gen.test$data

##  -------------------------
standardize<-function(mat){
  mat.std<-apply(mat,2,function(x) (x-mean(x))/sd(x))
  return(mat.std)
}


## simulation main body ----------------------
start<-proc.time()

library(doParallel) # foreach
cl<-makeCluster(10)
registerDoParallel(cl)

#---------------------------------------------------#
#  Configuration Panel #
#--------------------------------------------------#
K.all=3 #Set K = truth
lambda.all<-c(0.01,0.03,0.05,0.08,0.10,0.12,0.14,0.16,0.18,0.2) # penalty parameter
tune_pars=cbind(K=rep(K.all,each=length(lambda.all)),lambda=rep(lambda.all,length(K.all)))

robust="huberk" # robust method
alpha=0.3 # for elastic net
# 5 fold
library(caret)
set.seed(123)
flds <- createFolds(1:n, k = 5, list = TRUE, returnTrain = FALSE)
#---------------------------------------------------#

# consistency of estimates
par.est.all <- foreach(rep = 1:n.rep, .packages = c("mclust","OGClust")) %dopar% {
  # print(rep)
  df <- data.train[[rep]]
  X1 <- df$X1
  X2 <- df$X2
  X <- cbind(X1, X2)
  df[, 3:(NG + 2)] <- standardize(as.matrix(df[, 3:(NG + 2)])) # standardize
  G <- as.matrix(df[, 3:(NG + 2)])
  Y <- df$Y
  
  df_test <- data.test[[rep]]
  X1_test <- df_test$X1
  X2_test <- df_test$X2
  X_test <- cbind(X1_test, X2_test)
  df_test[, 3:(NG + 2)] <- standardize(as.matrix(df_test[, 3:(NG + 2)])) # standardize
  G_test <- as.matrix(df_test[, 3:(NG + 2)])
  Y_test <- df_test$Y
  
  onerep.res <- list()
  for (l in 1:nrow(tune_pars)) {
    # parameters
    K <- tune_pars[l, "K"]
    lambda <- tune_pars[l, "lambda"]
    # 2 initials and choose the best
    n.int <- 2
    est.par <- matrix(, ncol = ((K - 1) * (NG + 1) + np + K + 6), nrow = n.int)
    for (i in 1:n.int) {
      beta_int <- runif(np, 0, 3)
      gamma_int <- runif((K - 1) * (NG + 1), 0, 1)
      beta0_int <- runif(K, 0, 3)
      sigma2_int <- runif(1, 1, 3)
      theta_int <- c(beta_int, gamma_int, beta0_int, sigma2_int)
      fit.res <- fit.OGClust(
        n = n, K = K, np = np, NG = NG, lambda = lambda, alpha = alpha,
        G = G, Y = Y, X = X, theta_int = theta_int, robust = robust
      )
      est.par[i, ] <- fit.res$res
    }
    # choose the best par estimates which has the maximum ll
    best.par <- est.par[which.min(est.par[, (K - 1) * (NG + 1) + K + np + 5]), ]
    best.lambda <- best.par[(K - 1) * (NG + 1) + K + 8]
    beta1_est <- best.par[1]
    beta2_est <- best.par[2]
    gamma_est <- best.par[3:((K - 1) * (NG + 1) + 2)]
    miu_est <- best.par[((K - 1) * (NG + 1) + 3):((K - 1) * (NG + 1) + 2 + K)]
    sigma2_est <- best.par[((K - 1) * (NG + 1) + K + 3)]
    
    gamma_est_matrix <- matrix(gamma_est, ncol = K - 1, byrow = T)
    gamma_est_matrix <- cbind(gamma_est_matrix, 0)
    n_genes_selected=sum(apply(gamma_est_matrix[-1,],1,function(x) sum(x!=0)>0))
    n_true_genes_selected=sum(apply(gamma_est_matrix[2:16,],1,function(x) sum(x!=0)>0))
    #G <- cbind(1, G)
    pai_est <- sapply(1:K, function(k) {
      exp(cbind(1, G) %*% gamma_est_matrix[, k, drop = F]) / rowSums(exp(cbind(1, G) %*% gamma_est_matrix))
    })
    f_est <- sapply(1:K, function(x) (1 / sqrt(2 * pi * sigma2_est)) * exp(-(Y - miu_est[x] - X1 * beta1_est - X2 * beta2_est)^2 / (2 * sigma2_est)))
    
    # calculate the expected value of Z
    w_est <- sapply(1:K, function(k) {
      (pai_est[, k] * f_est[, k]) / diag(pai_est %*% t(f_est))
    })
    grp_assign_prd <- apply(w_est, 1, which.max)
    ARI <- adjustedRandIndex(grp_assign_prd, gen.train$grp_assign[rep, ])
    Y_prd <- miu_est[grp_assign_prd] + X1 * beta1_est + X2 * beta2_est
    
    # apply to validation set
    #validation group assignment and y prediction
    pai_est_test=sapply(1:K, function(k) 
      exp(cbind(1,G_test) %*% gamma_est_matrix[,k,drop=F])/(rowSums(exp(cbind(1,G_test) %*% gamma_est_matrix))))
    #Y_prd_test <- miu_est[grp_assign_prd] + X1 * beta1_est + X2 * beta2_est
    f_est<-sapply(1:K, function(x) (1/sqrt(2*pi*sigma2_est))*exp(-(Y_test-miu_est[x]-X1_test * beta1_est - X2_test * beta2_est)^2/(2*sigma2_est)))
    
    #calculate the expected value of Z
    w_est_test=sapply(1:K, function(k) 
      (pai_est_test[,k]*f_est[,k])/diag(pai_est_test %*% t(f_est)))
    grp_assign_prd_test<-apply(w_est_test,1,which.max)
    Y_prd_test <- miu_est[grp_assign_prd_test] + X1_test * beta1_est + X2_test * beta2_est
    ARI_test <- adjustedRandIndex(grp_assign_prd_test, gen.test$grp_assign[rep, ])
    onerep.res[[l]]<-list(K = K, best.lambda = best.lambda, best.par = best.par, prob = pai_est, post_prob = w_est, 
                     grp_assign_prd = grp_assign_prd, grp_assign_prd_test = grp_assign_prd_test,
                     ARI = ARI, ARI_test = ARI_test, Y_prd = Y_prd, Y_prd_test=Y_prd_test,
                     n_genes_selected=n_genes_selected, n_true_genes_selected=n_true_genes_selected)
  }
  return(onerep.res)
}

save(par.est.all, gen.train, gen.test, file = "ogclust_K3_1000genes_6groups_gamma1_delta3_lognormal_huberk.rdata")

ARI <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
ARI_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
n_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
n_true_markers<-matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2 <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
RMSE_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)
R2_test <- matrix(nrow=nrow(tune_pars),ncol=n.rep)

for(j in 1:nrow(tune_pars)){
  for (i in 1:n.rep) {
    #train
    Y <- gen.train$data[[i]]$Y
    Z <- gen.train$data[[i]]$Z
    X1 <- gen.train$data[[i]]$X1
    X2 <- gen.train$data[[i]]$X2
    Y_prd <- par.est.all[[i]][[j]]$Y_prd
    e=Y-(miu[Z]+beta1*X1+beta2*X2)
    RMSE[j,i]<- sqrt(sum((Y_prd - (Y-e))^2) / n)
    R2[j,i] <- 1 - sum((Y_prd-(Y-e))^2) / sum(((Y-e) - mean(Y))^2)
    ARI[j,i] <- par.est.all[[i]][[j]]$ARI
    n_markers[j,i]<- par.est.all[[i]][[j]]$n_genes_selected
    n_true_markers[j,i]<- par.est.all[[i]][[j]]$n_true_genes_selected    
    grp_assign_prd <- par.est.all[[i]][[j]]$grp_assign_prd
    # pvalues[i]=kruskal.test(Y ~ grp_assign_prd)$p.value
    print(table(grp_assign_prd,gen.train$grp_assign[i,]))
    
    #test
    Y_test <- gen.test$data[[i]]$Y
    Z_test <- gen.test$data[[i]]$Z
    X1_test <- gen.test$data[[i]]$X1
    X2_test <- gen.test$data[[i]]$X2
    Y_prd_test <- par.est.all[[i]][[j]]$Y_prd_test
    e_test=Y_test-(miu[Z_test]+beta1*X1_test+beta2*X2_test)
    RMSE_test[j,i]<- sqrt(sum((Y_prd_test - (Y_test-e_test))^2) / n)
    R2_test[j,i] <- 1 - sum((Y_prd_test-(Y_test-e_test))^2) / sum(((Y_test-e_test) - mean(Y_test-e_test))^2)
    ARI_test[j,i] <- par.est.all[[i]][[j]]$ARI_test
  }
}

ARI 
#[1] 0.8113398 0.8612458 0.8187639 0.8350720 0.7644884
#[6] 0.8353246 0.8441081 0.4216053 0.7630068 0.4806123
mean(ARI) #0.7435567
# summary(pvalues)
ARI_test

mean(ARI_test)

n_markers

n_true_markers

apply(RMSE,1,mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.072   2.188   2.265   2.409   2.312   3.235  
apply(R2,1,mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3327  0.5072  0.5291  0.5090  0.5641  0.6020 

apply(RMSE_test,1,mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.072   2.188   2.265   2.409   2.312   3.235  
apply(R2_test,1,mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3327  0.5072  0.5291  0.5090  0.5641  0.6020 