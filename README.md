# ogClust
An R package for outcome-guided clustering which identifies subgroups of outcome-associated patients towards precision medicine. 

## Installation
ogClust package files are in `package/` folder, You can install by copying and paste the following code into R

```
devtools::install_github("liupeng2117/ogClust/package")
```

Alternatively, download the `tar.gz` zipped file and install ogClust package using the code below

```
install.packages("~/ogClust_0.0.0.9000.tar.gz",repos=NULL,type="source")
```

## ogClust with continous outcome 
* Call the LGRC lung disease data. The data contains the lung gene expression data from 315 samples with 2000 genes, 3 prognostic covariates (`BMI`, `age` and `gender`) and one outcome `fev1pd`(FEV1 % predicted).
* `G` is a gene expression matrices with 315(samples) rows and 2000(genes) columns. `X` is a covariate matrix with 315 rows and 3(covariates) columns. `Y` is a numeric vector of a single continous outcome.

```r
library(ogClust)
data(lung)

G=lung$G # gene expression matrix
X=lung$X # covariates matrix
Y=lung$Y # outcome
n=nrow(G) # number of samples
NG=ncol(G) # number of genes
np=ncol(X) # number of covariates
```
   
* Set the number of clusters `K=3` and the tuning parameters(`lambda` and `alpha`)

```r
K=3
lambda=0.13
alpha=0.5
```

* Set the initial values for the parameters

```r
beta_int = runif(np, 0, 3)
gamma_int = runif((K - 1) * (NG + 1), 0, 1)
beta0_int = runif(K, 0, 3)
sigma2_int = runif(1, 1, 3)
theta_int = c(beta_int, gamma_int, beta0_int, sigma2_int)
``` 
    
* Fit ogClust model with continous outcome using `fit.ogClust` function
 
```r
fit.res<-fit.ogClust(n=n, K=K, np=np, NG=NG, lambda=lambda, 
                    alpha=alpha, G=G, Y=Y, X=X,theta_int=theta_int)
``` 

* Fitting result

`fit.res` is a list objects of class `ogClust` with four elements:
    + `res`: a vector of parameter estimates,likelihood `ll`, `R2`, `AIC`, `BIC` and tuning parameter `lambda`
    + `prob`: predicted probability for belonging to each subgroup
    + `Y_prd`: predicted outcome `Y`
    + `grp_assign`: group assignement

```r
str(fit.res)
#> List of 4
#>  $ res       : Named num [1:4014] 0.1917 0.0717 0.0113 0.8695 1.0684 ...
#>   ..- attr(*, "names")= chr [1:4014] "" "" "" "" ...
#>  $ prob      : num [1:315, 1:3] 0.0969 0.1019 0.2806 0.5142 0.8193 ...
#>  $ Y_prd     : num [1:315] 0.289 0.307 0.527 0.584 0.639 ...
#>  $ grp_assign: int [1:315] 3 2 2 1 1 2 3 3 3 3 ...
#>  - attr(*, "class")= chr "ogClust"
```
 
* Fit ogClust model with robust estimation, use `robust` option to choose which robust estimation procedure to use, default is `none`.
    + `robust = "median"`: use median-truncated loss in numerical estimation procedure 
    + `robust = "huber"`: use Huber loss in numerical estimation procedure 
    + `robust = "hubertf"`: use adaptive Huber method
 
```r
fit.res2<-fit.ogClust(n=n, K=K, np=np, NG=NG, lambda=lambda,
                       alpha=0.5, G=G, Y=Y, X=X,theta_int=theta_int,robust = "median")
fit.res3<-fit.ogClust(n=n, K=K, np=np, NG=NG, lambda=lambda,
                       alpha=0.5, G=G, Y=Y, X=X,theta_int=theta_int, robust ="huber",tau=4.345)
fit.res4<-fit.ogClust(n=n, K=K, np=np, NG=NG, lambda=lambda,
                       alpha=0.5, G=G, Y=Y, X=X,theta_int=theta_int, robust="hubertf")
```

## ogClust with survival outcome
* Call the example survival data. The example contains the gene expression data from 600 subjects with 1000 genes, 2 prognostic covariates (`X1` and `X2`) and one survival outcome.
* `G` is a gene expression matrices with 600(samples) rows and 2000(genes) columns. `X` is a covariate matrix with 600 rows and 2(covariates) columns. `Y` is a numeric vector of survival time, `delta` is a binary indicator of censoring.

```r
data("surv_dt")

G=surv_dt[,3:1002] # gene expression matrix
X=surv_dt[,1:2] # covariate matrix
Y=surv_dt$time # survival time 
delta=surv_dt$event # censoring indicator
n=nrow(G) # number of samples
NG=ncol(G) # number of genes
np=ncol(X) # number of covariates
```

* Set the number of clusters `K=3` and the tuning parameters(`lambda` and `alpha`)

```r
K=3
lambda=0.13
```

* Set the initial values for the parameters

```r
beta_int = runif(np, 0, 3)
gamma_int = runif((K - 1) * (NG + 1), 0, 1)
beta0_int = runif(K, 0, 3)
sigma2_int = runif(1, 1, 3)
theta_int = c(beta_int, gamma_int, beta0_int, sigma2_int)
```

* Fit ogClust model with survival outcome using `fit.ogClust.surv` function

```r
fit.res<-fit.ogClust.surv(n=n, K=K, np=np, NG=NG, lambda=lambda,
                            alpha=0.5, G=G, Y=Y, X=X, delta, theta_int=theta_int)
```
