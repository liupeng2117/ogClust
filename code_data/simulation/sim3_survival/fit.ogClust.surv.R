fit.OGClust.surv <- function(n, K, np, NG, lambda, alpha, G, Y, X, delta, theta_int) {
  G=cbind(1,G)
  theta_est = EM.surv(theta_int, lambda = lambda, n = n, G = G, Y = Y, X = X, delta=delta, np = np, K = K, NG = NG, alpha = alpha)$theta
  
  # estimated parameters
  beta_est = theta_est[1:np]
  gamma_est = theta_est[(np + 1):((K - 1) * (NG + 1) + np)]
  beta0_est = theta_est[((K - 1) * (NG + 1) + np + 1):((K - 1) * (NG + 1) + np + K)]
  sigma2_est = theta_est[((K - 1) * (NG + 1) + np + K + 1)]
  
  gamma_est_matrix = matrix(gamma_est, ncol = K - 1, byrow = T)
  gamma_est_matrix = cbind(gamma_est_matrix, 0)
  #G = cbind(1, G)
  pai_est = sapply(1:K, function(k) exp(G %*% gamma_est_matrix[, k, drop = F])/rowSums(exp(G %*% gamma_est_matrix)))
  #f_est <- sapply(1:K, function(x) (1/sqrt(2 * pi * sigma2_est)) * exp(-(Y - beta0_est[x] - X %*% beta_est)^2/(2 * sigma2_est)))
  f_est=sapply(1:K, function(x) f_calc(Y1=Y,X1=X,beta = beta_est,mu=beta0_est[x],sigma2 = sigma2_est,delta = delta))
  f_est=t(apply(f_est,1,function(x) x/sum(x)))
  idx=which.max(beta0_est)
  test=apply(f_est,1,sum)
  
  f_est[which(test<10^-3 | is.na(test)),idx]=1
  f_est[is.na(f_est)]=0
  
  (ll = sum(log(diag(pai_est %*% t(f_est)))))
  # calculate the expected value of Y and R2
  Y_prd = apply(sapply(1:K, function(x) pai_est[, x] * (beta0_est[x] + X %*% beta_est)), 1, sum)
  R2 = 1 - sum((Y - Y_prd)^2)/sum((Y - mean(Y))^2)
  
  # Calculate AIC BIC
  AIC = 2 * sum(theta_est != 0) - 2 * ll
  BIC = log(n) * sum(theta_est != 0) - 2 * ll
  #EBIC = BIC + 2 * (1 - 1/(2 * log(length(theta_est), base = n))) * log(choose(length(theta_est), sum(theta_est != 0)))
  
  # prosterior prob
  w_est = sapply(1:K, function(k) (pai_est[, k] * f_est[, k])/diag(pai_est %*% t(f_est)))
  cl.assign <- apply(w_est, 1, which.max)
  
  return(list(res = c(theta_est, ll, R2, AIC, BIC,lambda), theta_est = theta_est, Y_prd = Y_prd, grp_assign = cl.assign))
}


EM.surv<-function(theta, lambda, n, G, Y, X, delta, np, K, NG, alpha=0.5){
  #-----------------------------------------------#
  l=1
  repeat{
    #print(l)
    beta_old=theta[1:np]      # beta np times 1
    gamma_old=theta[(1+np):((K-1)*(NG+1)+np)]               ####?
    
    miu_old=theta[((K-1)*(NG+1)+np+1):((K-1)*(NG+1)+np+K)] ####?
    
    sigma2_old=theta[((K-1)*(NG+1)+np+K+1):length(theta)]
    
    
    #==E-STEP==#
    gamma_old_matrix=matrix(gamma_old,ncol=K-1,byrow=T)
    gamma_old_matrix=cbind(gamma_old_matrix,0)
    pai_old=sapply(1:K, function(k) 
      exp(G %*% gamma_old_matrix[,k,drop=F])/rowSums(exp(G %*% gamma_old_matrix)))
    #----------------------change-------------------------------------------------#
    
    
    
    
    #f_old<-sapply(1:K, function(x) (1/sqrt(2*pi*sigma2_old))*exp(-(Y-miu_old[x] - X %*% beta_old)^2/(2*sigma2_old)))#important
    
    f_old=sapply(1:K, function(x) f_calc(Y1=Y,X1=X,beta = beta_old,mu=miu_old[x],sigma2 = sigma2_old,delta = delta))
    f_old=t(apply(f_old,1,function(x) x/sum(x)))
    idx=which.max(miu_old)
    test=apply(f_old,1,sum)
    
    f_old[which(test<10^-3 | is.na(test)),idx]=1
    f_old[is.na(f_old)]=0
    #------------------------------------------------------------------------------#
    
    #calculate the expected value of Z
    w_old=sapply(1:K, function(k) 
      (pai_old[,k]*f_old[,k])/diag(pai_old %*% t(f_old)))
    #==M-STEP==#
    library(glmnet)
    gamma_new_matrix = tryCatch({
      fit<-glmnet(x=G[,-1], y=w_old, lambda=lambda, family = "multinomial", alpha = alpha, type.multinomial = "grouped")
      gamma_new_matrix=rbind(t(fit$a0),sapply(1:K, function(x) as.numeric(fit$beta[[x]])))
      gamma_new_matrix=sapply(1:K, function(x) gamma_new_matrix[,x]-gamma_new_matrix[,K])
    }, error =function(e) {
      return(gamma_old_matrix)
    })
    
    
    
    #---------------------------change----------------------------------#
    
    library(survival)
    dt0=data.frame(Y,delta,X)
    dt1=dt0
    for(k in 1:(K-1)) dt1 = rbind.data.frame(dt1,dt0)
    dt2=dt1
    for(k in 1:K) dt2=cbind.data.frame(dt2, rep(diag(K)[k,],each=dim(dt0)[1]))
    colnames(dt2)[(ncol(dt2)-K+1):ncol(dt2)]<-paste0("mu",1:K)
    weights=vector()
    for(k in 1:K){weights=c(weights,w_old[,k])}
    weights[which(weights==0)]=10^-3

      fit=survreg(eval(parse(text=paste("Surv(Y,delta)~-1+X1+X2+",paste(paste0("mu",1:K), collapse = " + ")))),
                  weights = weights,data = dt2,dist = "loglogistic", robust=TRUE) #
      miu_new=fit$coefficients[(np+1):(np+K)]
      sigma2_new=fit$scale
      #sigma2_new=1/150
      beta_new=fit$coefficients[1:np]
    
    theta_new=c(beta_new, 
                as.numeric(t(gamma_new_matrix[,-K])), 
                miu_new, sigma2_new)
    dis=sqrt(sum((theta_new-theta)^2,na.rm=T))
    theta=theta_new
    #print(theta)
    l=l+1
    if(dis<1*10^(-6) |l>500){
      break
    }
  }
  return(list(theta=theta,w=w_old,l=l))
}

#-------------------likelihood--------------------------#
f_calc=function(Y1,X1,beta,mu,sigma2,delta,K){
  Z=log(Y1)
  X0=X1%*%beta
  mu0=mu
  W=(Z-X0-mu0)/sigma2
  pdf_calc=function(delta,W,ind){
    #if(exp(W[ind])==Inf){return(1^(delta[ind]))}
    #else{
      return((1/sigma2)*(exp(W[ind])/(1+exp(W[ind]))^2)^(delta[ind]))
      #}
  }
  res=sapply(1:length(Y1),function(ind) (1/(1+exp(W[ind])))^(1-delta[ind])*pdf_calc(delta,W,ind))
  return(res)
}