#' Title Fit ogClust mixture model
#'
#' @param n the number of samples
#' @param K the number of subgroups
#' @param np the number of
#' @param NG the number of genes
#' @param lambda the penalty parameter for sparsity
#' @param alpha the tuning parameter for elastic net
#' @param G the input gene expression matrix
#' @param Y the input vector of outcome
#' @param X the input vector of covariates
#' @param theta_int initial values of parameters
#' @param robust choose the type of robustness
#' @param tau set the cutoff for huber loss if robust is huberk
#'
#' @details to be filled later
#'
#' @return An objects of class `ogClust`
#' \itemize{
#'  \item{'res'}{a vector of parameter estimates,likelihood `ll`, `R2`, `AIC`, `BIC` and tuning parameter `lambda`}
#'  \item{'prob'}{predicted probability for belonging to each subgroup}
#'  \item{'Y_prd'}{predicted outcome}
#'  \item{'grp_assign'}{prediced group assignement}
#' }
#' @export fit.ogClust
#' @importFrom glmnet glmnet
#' @importFrom tfHuber huberReg
#' @importFrom stats median
#'
#' @examples
#'   data(lung) #load lung dataset
#'
#'   # extract gene expression G, covariate X, survival time Y
#'   G=lung$G
#'   X=lung$X
#'   Y=lung$Y
#'
#'   # number of subjects
#'   n=nrow(G)
#'   # number of genes
#'   NG=ncol(G)
#'   # number of covariates
#'   np=ncol(X)
#'   # number of clusters
#'   K=3
#'   # tuning parameter
#'   lambda=0.13
#'
#'   # set initial values
#'   beta_int = runif(np, 0, 3)
#'   gamma_int = runif((K - 1) * (NG + 1), 0, 1)
#'   beta0_int = runif(K, 0, 3)
#'   sigma2_int = runif(1, 1, 3)
#'   theta_int = c(beta_int, gamma_int, beta0_int, sigma2_int)
#'
#'   # fit ogClust, robust is FALSE
#'   fit.res<-fit.ogClust(n=n, K=K, np=np, NG=NG, lambda=lambda,
#'                     alpha=0.5, G=G, Y=Y, X=X,theta_int=theta_int)
#'   # fit ogClust, robust method is median-truncation
#'   fit.res2<-fit.ogClust(n=n, K=K, np=np, NG=NG, lambda=lambda,
#'                      alpha=0.5, G=G, Y=Y, X=X,theta_int=theta_int,robust = 'median')
#'   # fit ogClust, robust method is Huber
#'   fit.res3<-fit.ogClust(n=n, K=K, np=np, NG=NG, lambda=lambda,
#'                      alpha=0.5, G=G, Y=Y, X=X,theta_int=theta_int, robust ='huber',tau=4.345)
#'   # fit ogClust, robust method is adaptive Huber
#'   fit.res4<-fit.ogClust(n=n, K=K, np=np, NG=NG, lambda=lambda,
#'                      alpha=0.5, G=G, Y=Y, X=X,theta_int=theta_int, robust='hubertf')
fit.ogClust <- function(n, K, np, NG, lambda, alpha, G, Y, X, theta_int, robust = "none", tau = 1.345) {
    stopifnot(robust %in% c("none", "huber", "median", "hubertf"))
    if (class(G) != "matrix") 
        G = as.matrix(G)
    if (class(X) != "matrix") 
        X = as.matrix(X)
    theta_est = EM(theta_int, lambda = lambda, n = n, G = G, Y = Y, X = X, np = np, K = K, NG = NG, alpha = alpha, robust = robust, tau = tau)
    
    # estimated parameters
    beta_est = theta_est[1:np]
    gamma_est = theta_est[(np + 1):((K - 1) * (NG + 1) + np)]
    beta0_est = theta_est[((K - 1) * (NG + 1) + np + 1):((K - 1) * (NG + 1) + np + K)]
    sigma2_est = theta_est[((K - 1) * (NG + 1) + np + K + 1)]
    
    gamma_est_matrix = matrix(gamma_est, ncol = K - 1, byrow = T)
    gamma_est_matrix = cbind(gamma_est_matrix, 0)
    G = cbind(1, G)
    pai_est = sapply(1:K, function(k) exp(G %*% gamma_est_matrix[, k, drop = F])/rowSums(exp(G %*% gamma_est_matrix)))
    f_est <- sapply(1:K, function(x) (1/sqrt(2 * pi * sigma2_est)) * exp(-(Y - beta0_est[x] - X %*% beta_est)^2/(2 * sigma2_est)))
    (ll = sum(log(diag(pai_est %*% t(f_est)))))
    # calculate the expected value of Y and R2
    Y_prd = apply(sapply(1:K, function(x) pai_est[, x] * (beta0_est[x] + X %*% beta_est)), 1, sum)
    R2 = 1 - sum((Y - Y_prd)^2)/sum((Y - mean(Y))^2)
    
    # Calculate AIC BIC
    AIC = 2 * sum(theta_est != 0) - 2 * ll
    BIC = log(n) * sum(theta_est != 0) - 2 * ll
    # EBIC = BIC + 2 * (1 - 1/(2 * log(length(theta_est), base = n))) * log(choose(length(theta_est), sum(theta_est != 0)))
    
    # prosterior prob
    w_est = sapply(1:K, function(k) (pai_est[, k] * f_est[, k])/diag(pai_est %*% t(f_est)))
    cl.assign <- apply(w_est, 1, which.max)
    final.res <- list(res = c(theta_est, ll = ll, R2 = R2, AIC = AIC, BIC = BIC, lambda = lambda), prob = pai_est, Y_prd = Y_prd, grp_assign = cl.assign)
    attr(final.res, "class") <- "ogClust"
    return(final.res)
}

EM <- function(theta, lambda, n, G, Y, X, np, K, NG, robust, alpha, tau = 1.345) {
    X = as.matrix(X)
    G = cbind(1, as.matrix(G))
    Y = Y
    
    l = 1
    repeat {
        # print(l)
        beta_old = theta[1:np]
        gamma_old = theta[(1 + np):((K - 1) * (NG + 1) + np)]
        miu_old = theta[((K - 1) * (NG + 1) + np + 1):((K - 1) * (NG + 1) + np + K)]
        sigma2_old = theta[((K - 1) * (NG + 1) + np + K + 1)]
        
        # ==E-STEP==#
        gamma_old_matrix = matrix(gamma_old, ncol = K - 1, byrow = T)
        gamma_old_matrix = cbind(gamma_old_matrix, 0)
        pai_old = sapply(1:K, function(k) exp(G %*% gamma_old_matrix[, k, drop = F])/rowSums(exp(G %*% gamma_old_matrix)))
        f_old <- sapply(1:K, function(x) (1/sqrt(2 * pi * sigma2_old)) * exp(-(Y - miu_old[x] - X %*% beta_old)^2/(2 * sigma2_old)))
        
        
        # calculate the expected value of Z
        w_old = sapply(1:K, function(k) (pai_old[, k] * f_old[, k])/diag(pai_old %*% t(f_old)))
        # ==M-STEP==#
        gamma_new_matrix = tryCatch({
            fit <- glmnet::glmnet(x = G[, -1], y = w_old, lambda = lambda, family = "multinomial", alpha = alpha, type.multinomial = "grouped")
            gamma_new_matrix = rbind(t(fit$a0), sapply(1:K, function(x) as.numeric(fit$beta[[x]])))
            gamma_new_matrix = sapply(1:K, function(x) gamma_new_matrix[, x] - gamma_new_matrix[, K])
        }, error = function(e) {
            return(gamma_old_matrix)
        })
        
        #---- non robust ----#
        if (robust == "none") {
            # update miu, alpha
            miu_new = sapply(1:K, function(k) {
                sum(w_old[, k] * (Y - X %*% beta_old), na.rm = T)/sum(w_old[, k], na.rm = T)
            })
            # update beta
            beta_new <- beta_old
            for (i in 1:np) {
                beta_new[i] <- sum(sapply(1:K, function(k) w_old[, k] * X[, i] * (Y - miu_new[k] - X[, -i, drop = F] %*% beta_new[-i])), na.rm = T)/sum(w_old * 
                  X[, i]^2, na.rm = T)
            }
            # update sigma2
            sigma2_new = sum(sapply(1:K, function(k) w_old[, k] * (Y - miu_new[k] - X %*% beta_new)^2), na.rm = T)/n
        }
        
        #--- median truncated ---#
        if (robust == "median") {
            e = sapply(1:K, function(k) Y - miu_old[k] - X %*% beta_old)
            # update miu, alpha
            miu_new = sapply(1:K, function(k) sum((w_old[, k] * (Y - X %*% beta_old))[abs(e[, k]) <= median(abs(e[, k]))], na.rm = T)/sum(w_old[abs(e[, k]) <= 
                median(abs(e[, k])), k], na.rm = T))
            # update beta1 beta2
            beta_new <- beta_old
            for (i in 1:np) {
                beta_new[i] = sum(unlist(sapply(1:K, function(k) (w_old[, k] * X[, i] * (Y - miu_new[k] - X[, -i, drop = F] %*% beta_new[-i]))[abs(e[, k]) <= 
                  median(abs(e[, k]))])), na.rm = T)/sum(unlist(sapply(1:K, function(k) (w_old[, k] * X[, i]^2)[abs(e[, k]) <= median(abs(e[, k]))])), na.rm = T)
            }
            # update sigma2
            sigma2_new = sum(unlist(sapply(1:K, function(k) (w_old[, k] * (Y - miu_new[k] - X %*% beta_new)^2)[abs(e[, k]) <= median(abs(e[, k]))])), na.rm = T)/sum(unlist(sapply(1:K, 
                function(k) w_old[abs(e[, k]) <= median(abs(e[, k])), k])), na.rm = T)
            
        }
        
        #--- huber ---#
        if (robust == "huber") {
            e = sapply(1:K, function(k) Y - miu_old[k] - X %*% beta_old)
            # update miu, alpha
            miu_new = sapply(1:K, function(k) sum(c((w_old[, k] * (Y - X %*% beta_old))[abs(e[, k]) <= tau], (w_old[, k] * tau * sign(e[, k]))[abs(e[, k]) > 
                tau]), na.rm = T)/sum(w_old[abs(e[, k]) <= tau, k], na.rm = T))
            # update beta
            beta_new = beta_old
            for (i in 1:np) {
                beta_new[i] = sum(sapply(1:K, function(k) c((w_old[, k] * X[, i] * (Y - miu_new[k] - X[, -i, drop = F] %*% beta_new[-i]))[abs(e[, k]) <= tau], 
                  (w_old[, k] * tau * sign(e[, k]))[abs(e[, k]) > tau])), na.rm = T)/sum(sapply(1:K, function(k) sum((w_old[, k] * X[, i]^2)[abs(e[, k]) <= tau], 
                  na.rm = T)), na.rm = T)
            }
            # update sigma2
            sigma2_new = sum(sapply(1:K, function(k) c((w_old[, k] * (Y - miu_new[k] - X %*% beta_new)^2)[abs(e[, k]) <= tau], (w_old[, k] * (2 * tau * abs(Y - 
                miu_new[k] - X %*% beta_new) - tau^2))[abs(e[, k]) > tau])), na.rm = T)/sum(sapply(1:K, function(k) sum(w_old[abs(e[, k]) <= tau, ], na.rm = T)), 
                na.rm = T)
        }
        
        #--- adaptive huber ---#
        if (robust == "hubertf") {
            beta_new = beta_old
            miu_new = miu_old
            grp.id <- t(1 * apply(w_old, 1, function(x) x == max(x)))[, -1]
            X_tf <- cbind(grp.id, X)
            
            listHuber = tfHuber::huberReg(X_tf, Y)
            mm <- diag(K)
            mm[1, ] <- 1
            miu_new <- listHuber$theta[1:K] %*% mm
            beta_new <- listHuber$theta[-c(1:K)]
            tau = listHuber$tauCoef
            e = sapply(1:K, function(k) Y - miu_new[k] - X %*% beta_new)
            sigma2_new = sum(sapply(1:K, function(k) c((w_old[, k] * (Y - miu_new[k] - X %*% beta_new)^2)[abs(e[, k]) <= tau], (w_old[, k] * (2 * tau * abs(Y - 
                miu_new[k] - X %*% beta_new) - tau^2))[abs(e[, k]) > tau])), na.rm = T)/sum(sapply(1:K, function(k) sum(w_old[abs(e[, k]) <= tau, ], na.rm = T)), 
                na.rm = T)
        }
        
        theta_new = c(beta_new, as.numeric(t(gamma_new_matrix[, -K])), miu_new, sigma2_new)
        dis = sqrt(sum((theta_new - theta)^2, na.rm = T))
        theta = theta_new
        l = l + 1
        if (dis < 1 * 10^(-6) | l > 500) {
            break
        }
    }
    return(theta)
}

