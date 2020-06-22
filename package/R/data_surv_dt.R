#' Simulated survival data
#'
#' An simulated survival data, with expression of 1000
#' genes, two prognostic covariates, survival time,
#' plus an indicator for censoring.
#'
#' @docType data
#'
#' @usage data(surv_dt)
#'
#' @format a list containing gene expression matrix `G`,
#' covariate matrix `X` and survival time `Y` and a binary
#' indicator `delta` for censoring.
#'
#' @author Peng Liu, 2020-06-21
#'
#' @keywords datasets
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
"surv_dt"
