## test function in test1.R
test_that("test fit.ogClust.surv function", {
  data("surv_dt")
  G=surv_dt[,3:1002]
  X=surv_dt[,1:2]
  Y=surv_dt$time
  delta=surv_dt$event
  n=nrow(G)
  NG=ncol(G)
  np=ncol(X)
  K=3
  lambda=0.13
  beta_int = runif(np, 0, 3)
  gamma_int = runif((K - 1) * (NG + 1), 0, 1)
  beta0_int = runif(K, 0, 3)
  sigma2_int = runif(1, 1, 3)

  theta_int = c(beta_int, gamma_int, beta0_int, sigma2_int)
  fit.res<-fit.ogClust.surv(n=n, K=K, np=np, NG=NG, lambda=lambda,
                            alpha=0.5, G=G, Y=Y, X=X, delta, theta_int=theta_int)
  expect_equal(length(fit.res), 4)
})

