rm(list = ls())
if (!require("MASS")) install.packages("MASS")
set.seed(100)


##### DGP #####
DGP <- function(K, T_, N, range_x, beta_i, alpha_i, mu_u_i, sigma_u_i){
  if(length(beta_i) != K){
    warning("length of beta_i != K")
    return(c())
  }
  if(length(mu_u_i) != T_){
    warning("length of mu_u_i != T")
    return(c())
  }
  if(length(diag(sigma_u_i)) != T_){
    warning("dims of sigma_u_i != (T,T)")
    return(c())
  }
  X <- lapply(as.list(1:N), function(i) array(runif(K*T_,range_x[1],range_x[2]), dim = c(T_, K)))
  U <- lapply(as.list(1:N), function(i) mvrnorm(n = 1, mu = mu_u_i, Sigma = sigma_u_i))
  Y <- mapply(function(X_i, u_i) alpha_i + X_i %*% beta_i + u_i, X, U, SIMPLIFY=F)
  # X <- list()
  # Y <- list()
  # for(i in 1:N){
  #     X_i <- array(runif(K*T_,range_x[1],range_x[2]), dim = c(T_, K))
  #     u_i = mvrnorm(n = 1, mu = mu_u_i, Sigma = sigma_u_i)
  #     Y_i = alpha_i + X_i %*% beta_i + u_i
  #     
  #     X[[i]] = X_i
  #     Y[[i]] = Y_i
  # }
  return(list(X = X, Y = Y))
}

OLS_FE <- function(X, Y){
  if(length(dim(X[[1]])) != 2){
    warning("X_i is not 2d array")
    return(c())
  }
  
  N = length(X)
  K = dim(X[[1]])[2]
  T_ = dim(X[[1]])[1]
  A = array(0, dim = c(K,K))
  B = rep(0, K)
  for(i in 1:N){
    x_i_mean = colMeans(X[[i]])
    y_i_mean = mean(Y[[i]])
    for(t in 1:T_){
      A = A + (X[[i]][t,] - x_i_mean) %*% t(X[[i]][t,] - x_i_mean)
      B = B + (X[[i]][t,] - x_i_mean) * (Y[[i]][t] - y_i_mean)
    }
  }
  beta_hat_fe = solve(A) %*% B
  return(beta_hat_fe)
}

dgp = DGP(K = 3,
          T_ = 3,
          N = 400,
          range_x = c(0,1),
          beta_i = c(1,2,3),
          alpha_i = 1,
          mu_u_i = c(0,0,0),
          sigma_u_i = diag(c(1,1,1)))
OLS_FE(dgp$X, dgp$Y)

