rm(list=ls())
source("~/subtyping/simulation/code/gen.data.K3.R")
library(OGClust)
# Simulation ------------------------------------------------
## Data generation --------------------
#============================
# 600 subjects, 2 subgroups, 6 clusters, 2 clinical variables x1, x2
# 100 genes G1 to G100
# G1 and G10 form two clusters and related to outcome
# G11 and G25 forms another three clusters but not related
# G26 to G100 are noise
# 10% outliers on Y
#============================
set.seed(123)
n.rep=100
NG=1000
n=600 # sample
np=2 # 2 covariates
beta1=1 # coefficient of clinical var1
beta2=1 # coefficient of clinical var2
gamma1=c(1,1,1,1,1,0,0,0,0,0,-1,-1,-1,-1,-1) # 1st set coefficient of genes 1 to 15, others are 0
gamma2=c(0,0,0,0,0,1,1,1,1,1,-1,-1,-1,-1,-1) # 2nd set coefficient of genes 1 to 15, others are 0
miu=c(1,6,11) # overall mean
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

gen.train<-gen.data.K3(n.rep, NG, n, theta, cluster_index, exp=0.5)
data.train<-gen.train$data

gen.test<-gen.data.K3(n.rep, NG, n, theta, cluster_index, exp=0.5)
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
K.all=c(2,3,6,9) #Set K = truth
lambda.all<-c(0.02,0.03,0.05,0.1,0.15) # penalty parameter
tune_pars=cbind(K=rep(K.all,each=length(lambda.all)),lambda=rep(lambda.all,length(K.all)))

robust="none" # robust method
alpha=0.3 # for elastic net
# 5 fold
#library(caret)
#set.seed(123)
#flds <- createFolds(1:n, k = 5, list = TRUE, returnTrain = FALSE)
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
  
  res.all <- list()
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
    res.all[[l]] <- best.par
  }
  
  BIC.all <- c()
  for (l in 1:nrow(tune_pars)) {
    K <- tune_pars[l, "K"]
    BIC.all[l] <- res.all[[l]][(K - 1) * (NG + 1) + np + K + 5]
  }
  
  print("All K and lambda DONE")
  
  # best K
  K <- best.K <- tune_pars[which.min(BIC.all), "K"]
  best.lambda <- tune_pars[which.min(BIC.all), "lambda"]
  print(paste0("best K=", best.K))
  print(paste0("best lambda=", best.lambda))
  
  best.par <- res.all[[which.min(BIC.all)]]
  print(best.par[c(1:2, 104:106)])
  best.lambda <- best.par[(K - 1) * (NG + 1) + K + 8]
  beta1_est <- best.par[1]
  beta2_est <- best.par[2]
  beta_est=c(beta1_est,beta2_est)
  gamma_est <- best.par[3:((K - 1) * (NG + 1) + 2)]
  miu_est <- best.par[((K - 1) * (NG + 1) + 3):((K - 1) * (NG + 1) + 2 + K)]
  sigma2_est <- best.par[((K - 1) * (NG + 1) + K + 3)]
  
  gamma_est_matrix <- matrix(gamma_est, ncol = K - 1, byrow = T)
  gamma_est_matrix <- cbind(gamma_est_matrix, 0)
  
  G <- cbind(1, G)
  pai_est <- sapply(1:K, function(k) {
    exp(G %*% gamma_est_matrix[, k, drop = F]) / rowSums(exp(G %*% gamma_est_matrix))
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
  Y_prd_test <- miu_est[grp_assign_prd] + X1 * beta1_est + X2 * beta2_est
  f_est<-sapply(1:K, function(x) (1/sqrt(2*pi*sigma2_est))*exp(-(Y_test-miu_est[x]-X_test %*% beta_est)^2/(2*sigma2_est)))
  
  #calculate the expected value of Z
  w_est_test=sapply(1:K, function(k) 
    (pai_est_test[,k]*f_est[,k])/diag(pai_est_test %*% t(f_est)))
  grp_assign_prd_test<-apply(w_est_test,1,which.max)
  ARI_test <- adjustedRandIndex(grp_assign_prd_test, gen.test$grp_assign[rep, ])
  
  
  return(list(K = K, best.lambda = best.lambda, best.par = best.par, prob = pai_est, post_prob = w_est,
              prob_test = pai_est_test, post_prob_test = w_est_test,
              grp_assign_prd = grp_assign_prd,grp_assign_prd_test = grp_assign_prd_test,
              ARI = ARI,ARI_test = ARI_test, Y_prd = Y_prd, Y_prd_test = Y_prd_test))
}

save(par.est.all, gen.train, gen.test, file = "ogclust_K3_1000genes_9groups_exp0.5_gamma1_delta5.rdata")

ARI <- c()
# pvalues<-c()
for (i in 1:n.rep) {
  ARI[i] <- par.est.all[[i]]$ARI
  Y <- gen.train$data[[i]]$Y
  grp_assign_prd <- par.est.all[[i]]$grp_assign_prd
  # pvalues[i]=kruskal.test(Y ~ grp_assign_prd)$p.value
  print(table(grp_assign_prd,gen.train$grp_assign[i,]))
}
ARI 
#[1] 0.9897129 0.9802020 0.9899096 0.9752050 0.9950593 0.9907367 0.9655690 0.9802408 0.9900802 0.9950613
median(ARI) #[1] 0.9898113

est.K<-c()
for (i in 1:n.rep) {
  est.K[i] <- par.est.all[[i]]$K
}
table(est.K)
#est.K
#3 
#10 

# Genes: FNs FPs
TPs<-c()
TNs<-c()
FNs<-c()
FPs<-c()
for(i in 1:100){
  best.par<-par.est.all[[i]]$best.par
  K<-est.K[i]
  gamma_est <- best.par[3:((K - 1) * (NG + 1) + 2)]
  gamma_est_matrix <- matrix(gamma_est, ncol = K - 1, byrow = T)
  gamma_est_matrix <- cbind(gamma_est_matrix, 0)
  gamma_est_matrix <- gamma_est_matrix[-1,]
  TPs[i]<-sum(which(apply(gamma_est_matrix, 1, function(x) sum(abs(x)>0)>=1)) %in% 1:15)
  TNs[i]<-sum(which(apply(gamma_est_matrix, 1, function(x) sum(abs(x)>0)==0)) %in% 16:1000)
  FNs[i]<-sum(which(apply(gamma_est_matrix, 1, function(x) sum(abs(x)>0)==0)) %in% 1:15)
  FPs[i]<-sum(which(apply(gamma_est_matrix, 1, function(x) sum(abs(x)>0)>=1)) %in% 16:1000)
}

summary(FNs)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0       0       0       0       0       0 
summary(FPs)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0       0       0       0       0       0 

RMSE <- c()
R2 <- c()
for (rep in 1:n.rep) {
  Y <- gen.train$data[[rep]]$Y
  Y_prd <- par.est.all[[rep]]$Y_prd
  RMSE[rep] <- sqrt(sum((Y_prd - Y)^2) / n)
  R2[rep] <- 1 - sum((Y - Y_prd)^2) / sum((Y - mean(Y))^2)
}
summary(RMSE)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.603   1.642   1.772   1.744   1.821   1.883 
summary(R2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5319  0.5397  0.5713  0.5863  0.6365  0.6559 