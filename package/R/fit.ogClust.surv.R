#' Title Fit ogClust mixture model
#'
#' @param n an integer defines number of samples
#' @param K an integer defines the number of subgroups
#' @param np an integer defines the number of covariates
#' @param NG an integer defines the number of features of omics data
#' @param lambda the regularization tuning parameter for sparsity
#' @param alpha the L2 regularization tuning parameter
#' @param G the matrix for omics data. The rows are samples and columns are features.
#' @param X the vector of covariates
#' @param Y the vector of time variable
#' @param delta a binary indicator of censoring for time-to-event outcome. 0 means censored and 1 means event observed.
#' @param theta_int a vector of initial values for parameters to estimate
#' @param dist distribution of the survival time, defualt is \code{"loglogistic"}. The choices are \code{"weibull"},
#' \code{"exponential"}, \code{"gaussian"}, \code{"logistic"},\code{"lognormal"} and \code{"loglogistic"}.
#'
#' @details The ogClust is a unified latent generative model to perform clustering constructed
#' from omics data \code{G} with the guidance of outcome \code{Y}, and with covariate \code{X} to account for
#' the variability that is not related to subgrouping. A modified EM algorithm is applied for
#' numeric computation such that the liklihood is maximized. A posterior probability is obtain
#' for each subject belonging to each cluster.
#'
#' ogClust method performs feature selection, latent subtype characterization and outcome prediction simultaneously.
#' We use either LASSO penalty or LASSO penalty plus L2 regularization( \eqn{lambda*(L1 + alpha*L2)})
#' Parameter \code{lambda} is the penalty tuning parameter, and \code{alpha} tunes the L2 regularization.
#' Time variable \code{Y} and censoring indicator \code{delta} together defines a single time-to-event outcome.
#' We use accelerated-failure time(AFT) model to model time-to-event outcome, \code{dist} defines
#' the distribution of survival time.
#'
#' @return An object with class \code{"ogClust"}
#' \itemize{
#'  \item{\code{par}}{ a list of parameter estimates}
#'  \item{\code{ll}}{likelihood}
#'  \item{\code{AIC}}{AIC}
#'  \item{\code{BIC}}{BIC}
#'  \item{\code{lambda}}{lambda}
#'  \item{\code{prob}}{ predicted probability for belonging to each subgroup}
#'  \item{\code{Y_prd}}{ predicted outcome}
#'  \item{\code{grp_assign}}{ prediced group assignement}
#' }
#' @export fit.ogClust.surv
#' @importFrom glmnet glmnet
#' @import survival
#'
#' @examples
#'   data('surv_dt') #load simulated survial data
#'
#'   # extract gene expression G, covariate X, survival time Y
#'   # and censoring indicator delta
#'   G=surv_dt[,3:1002]
#'   X=surv_dt[,1:2]
#'   Y=surv_dt$time
#'   delta=surv_dt$event
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
#'   # fit ogClust
#'   fit.res<-fit.ogClust.surv(n=n, K=K, np=np, NG=NG, lambda=lambda,
#'                          alpha=0.5, G=G, Y=Y, X=X, delta, theta_int=theta_int)
fit.ogClust.surv <- function(n, K, np, NG, lambda, alpha, G, Y, X, delta, theta_int, dist = "loglogistic") {
    G = cbind(1, as.matrix(G))
    X = as.matrix(X)
    theta_est = EM.surv(theta_int, lambda = lambda, n = n, G = G, Y = Y, X = X, delta = delta, np = np, K = K, NG = NG, alpha = alpha, dist = dist)$theta

    # estimated parameters
    beta_est = theta_est[1:np]
    gamma_est = theta_est[(np + 1):((K - 1) * (NG + 1) + np)]
    beta0_est = theta_est[((K - 1) * (NG + 1) + np + 1):((K - 1) * (NG + 1) + np + K)]
    sigma2_est = theta_est[((K - 1) * (NG + 1) + np + K + 1)]

    gamma_est_matrix = matrix(gamma_est, ncol = K - 1, byrow = T)
    gamma_est_matrix = cbind(gamma_est_matrix, 0)
    
    par_est<-list(beta0=beta0_est, beta=beta_est, sigma2=sigma2_est,gamma=gamma_est_matrix)
    pai_est = sapply(1:K, function(k) exp(G %*% gamma_est_matrix[, k, drop = F])/rowSums(exp(G %*% gamma_est_matrix)))
    f_est = sapply(1:K, function(x) f_calc(Y1 = Y, X1 = X, beta = beta_est, mu = beta0_est[x], sigma2 = sigma2_est, delta = delta))
    f_est = t(apply(f_est, 1, function(x) x/sum(x)))
    idx = which.max(beta0_est)
    test = apply(f_est, 1, sum)

    f_est[which(test < 10^-3 | is.na(test)), idx] = 1
    f_est[is.na(f_est)] = 0

    (ll = sum(log(diag(pai_est %*% t(f_est)))))

    # Calculate AIC BIC
    AIC = 2 * sum(theta_est != 0) - 2 * ll
    BIC = log(n) * sum(theta_est != 0) - 2 * ll

    # prosterior prob
    w_est = sapply(1:K, function(k) (pai_est[, k] * f_est[, k])/diag(pai_est %*% t(f_est)))
    cl.assign <- apply(w_est, 1, which.max)
    
    # calculate the expected value of Y and R2
    #Y_prd = apply(sapply(1:K, function(x) pai_est[, x] * (beta0_est[x] + X %*% beta_est)), 1, sum) #soft assignment
    Y_prd = beta0_est[cl.assign] + X %*% beta_est #hard assignment
    #R2 = 1 - sum((Y - Y_prd)^2)/sum((Y - mean(Y))^2)
    
    final.res <- list(par=par_est, ll = ll, AIC = AIC, BIC = BIC, lambda = lambda,
                      Y_prd=Y_prd, grp_assign=cl.assign)
    attr(final.res, "class") <- "ogClust"
    return(final.res)
}


EM.surv <- function(theta, lambda, n, G, Y, X, delta, np, K, NG, alpha = 0.5, dist) {
    #-----------------------------------------------#
    l = 1
    repeat {
        beta_old = theta[1:np]
        gamma_old = theta[(1 + np):((K - 1) * (NG + 1) + np)]

        miu_old = theta[((K - 1) * (NG + 1) + np + 1):((K - 1) * (NG + 1) + np + K)]

        sigma2_old = theta[((K - 1) * (NG + 1) + np + K + 1):length(theta)]


        # ==E-STEP==#
        gamma_old_matrix = matrix(gamma_old, ncol = K - 1, byrow = T)
        gamma_old_matrix = cbind(gamma_old_matrix, 0)
        #pai_old = sapply(1:K, function(k) exp(G %*% gamma_old_matrix[, k, drop = F])/rowSums(exp(G %*% gamma_old_matrix)))
        pai_old = sapply(1:K, function(k) 1/rowSums(exp(sweep(G %*% gamma_old_matrix, 1, G %*% gamma_old_matrix[, k, drop = F]))))

        f_old = sapply(1:K, function(x) f_calc(Y1 = Y, X1 = X, beta = beta_old, mu = miu_old[x], sigma2 = sigma2_old, delta = delta))
        f_old = t(apply(f_old, 1, function(x) x/sum(x)))
        idx = which.max(miu_old)
        test = apply(f_old, 1, sum)

        f_old[which(test < 10^-3 | is.na(test)), idx] = 1
        f_old[is.na(f_old)] = 0
        # calculate the expected value of Z
        w_old = sapply(1:K, function(k) (pai_old[, k] * f_old[, k])/diag(pai_old %*% t(f_old)))

        # ==M-STEP==#
        #gamma_new_matrix = tryCatch({

        if(any(apply(w_old,2, function(x) sum(abs(x)))<1e-05)){
            w_old[,apply(w_old,2, function(x) sum(abs(x)))<1e-05]<-1e-05/n
        }
            fit <- glmnet::glmnet(x = G[, -1], y = w_old, lambda = lambda, family = "multinomial", alpha = alpha, type.multinomial = "grouped")
            gamma_new_matrix = rbind(t(fit$a0), sapply(1:K, function(x) as.numeric(fit$beta[[x]])))
            gamma_new_matrix = sapply(1:K, function(x) gamma_new_matrix[, x] - gamma_new_matrix[, K])
        #}, error = function(e) {
        #    return(gamma_old_matrix)
        #})

        if (is.null(colnames(X))) {
            colnames(X) = paste0("X", 1:np)
        }
        dt0 = data.frame(Y, delta, X)
        dt1 = dt0
        for (k in 1:(K - 1)) dt1 = rbind.data.frame(dt1, dt0)
        dt2 = dt1
        for (k in 1:K) dt2 = cbind.data.frame(dt2, rep(diag(K)[k, ], each = dim(dt0)[1]))
        colnames(dt2)[(ncol(dt2) - K + 1):ncol(dt2)] <- paste0("mu", 1:K)
        weights = vector()
        for (k in 1:K) {
            weights = c(weights, w_old[, k])
        }
        weights[which(weights == 0)] = 10^-3

        fit = survival::survreg(eval(parse(text = paste("Surv(Y,delta)~-1", paste(colnames(X), collapse = " + "), paste(paste0("mu", 1:K), collapse = " + "),
            sep = " + "))), weights = weights, data = dt2, dist = dist, robust = FALSE)
        miu_new = fit$coefficients[(np + 1):(np + K)]
        sigma2_new = fit$scale
        beta_new = fit$coefficients[1:np]

        theta_new = c(beta_new, as.numeric(t(gamma_new_matrix[, -K])), miu_new, sigma2_new)
        dis = sqrt(sum((theta_new - theta)^2, na.rm = T))
        theta = theta_new
        l = l + 1
        if (dis < 1 * 10^(-6) | l > 500) {
            break
        }
    }
    return(list(theta = theta, w = w_old, l = l))
}

#-------------------likelihood--------------------------#
f_calc = function(Y1, X1, beta, mu, sigma2, delta, K) {
    Z = log(Y1)
    X0 = X1 %*% beta
    mu0 = mu
    W = (Z - X0 - mu0)/sigma2
    pdf_calc = function(delta, W, ind) {
        return((1/sigma2) * (exp(W[ind])/(1 + exp(W[ind]))^2)^(delta[ind]))
    }
    res = sapply(1:length(Y1), function(ind) (1/(1 + exp(W[ind])))^(1 - delta[ind]) * pdf_calc(delta, W, ind))
    return(res)
}
