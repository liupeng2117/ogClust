#' LGRC lung disease data
#'
#' Gene expression data are collected from Gene Expression Omnibus (GEO) GSE47460
#' and clinical information obtained from Lung Genomics Research Consortium
#' (https://ltrcpublic.com/). The data has n=315 pateints diagnosed by
#' two most representative lung disease subtypes: chronic obstructive pulmonary disease (COPD)
#' and interstitial lung disease (ILD), expression of 2000 genes, and values of
#' three prognostic covariates.
#' The outcome is FEV1%prd, which is a person's measured FEV1(the
#' volume of air a person can exhale during the first second of
#' forced expiration) normalized by the predicted FEV1 with healthy lung.
#'
#' @docType data
#'
#' @usage data(lung)
#'
#' @format a list containing gene expression matrix `G`, covariate matrix `X` and outcome `Y`
#'
#' @keywords datasets
#'
#' @source Lung Genomics Research Consortium (https://ltrcpublic.com/)
#'         Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/)
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
"lung"
