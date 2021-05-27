test_that("test fit.ogClust function", {
  data(lung)
  G=lung$G
  X=lung$X
  Y=lung$Y
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
  fit.res<-fit.ogClust(n=n, K=K, np=np, NG=NG, lambda=lambda,
                       alpha=0.5, G=G, Y=Y, X=X,theta_int=theta_int)
  fit.res2<-fit.ogClust(n=n, K=K, np=np, NG=NG, lambda=lambda,
                       alpha=0.5, G=G, Y=Y, X=X,theta_int=theta_int,robust = "median")
  fit.res3<-fit.ogClust(n=n, K=K, np=np, NG=NG, lambda=lambda,
                       alpha=0.5, G=G, Y=Y, X=X,theta_int=theta_int, robust ="huber",tau=4.345)
  fit.res4<-fit.ogClust(n=n, K=K, np=np, NG=NG, lambda=lambda,
                       alpha=0.5, G=G, Y=Y, X=X,theta_int=theta_int, robust="hubertf")
  expect_equal(length(fit.res), 4)
  expect_equal(length(fit.res2), 4)
  expect_equal(length(fit.res3), 4)
  expect_equal(length(fit.res4), 4)
})
