rm(list = ls())
if (!require("MASS")) install.packages("MASS")
if (!require("ggplot2")) install.packages("ggplot2")
library(purrr)
library(cowplot)
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
  return(list(beta_hat = beta_hat_fe))
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


##### simulation #####


# 
sim <- function(K, beta_true, all_n, all_T, nsims){
  beta_hat_1 = array(NA, dim = c(length(all_n), length(all_T), nsims))
  mean_beta_hat_1 = array(NA, dim = c(length(all_n), length(all_T)))
  var_beta_hat_1 = array(NA, dim = c(length(all_n), length(all_T)))
  for(i in 1:length(all_n)){
    N = all_n[i]
    for(j in 1:length(all_T)){
      T_ = all_T[j]
      for(h in 1:nsims){
        dgp = DGP(K = K,
                  T_ =T_,
                  N = N,
                  range_x = c(0,1),
                  beta_i = beta_true,
                  alpha_i = 1,
                  mu_u_i = rep(0, T_),
                  sigma_u_i = diag(rep(1, T_)))
        result = OLS_FE(dgp$X, dgp$Y)
        beta_hat = result$beta_hat
        beta_hat_1[i, j, h] = beta_hat[1]
      }
      mean_beta_hat_1[i,j] = mean(beta_hat_1[i, j,])
      var_beta_hat_1[i,j] = var(beta_hat_1[i, j,])
    }
  }
  return(list(beta_hat_1=beta_hat_1, mean_beta_hat_1=mean_beta_hat_1, var_beta_hat_1=var_beta_hat_1))
}

##### sim1 #####
sim1 <- sim(K = 3, beta_true = c(1:3), all_n = c(1000, 2000, 3000), all_T = c(3), nsims = 100)
beta_hat_1 = sim1$beta_hat_1
beta_true_1 = 1

# generate point and box plot for every data size
point_plot.list <- list()
box_plot.list <- list()
for(i in 1:length(all_n)){
  N <- all_n[i]
  plot.heigth <- max(beta_true_1 - min(beta_hat_1), max(beta_hat_1) - beta_true_1)
  point_plot.list[[i]] <- ggplot() +
    geom_point(aes(y = beta_hat_1[i,1,], x = c(1:nsims))) +
    geom_hline(yintercept = beta_true_1, color = I("black")) +
    ylim(beta_true_1 - plot.heigth, beta_true_1 + plot.heigth) +
    labs(subtitle = paste0("sample size = ", N), x = "iteration", y = "beta_hat_1")
  
  box_plot.list[[i]] <- ggplot() +
    geom_boxplot(aes(y = beta_hat_1[i,1,])) +
    ylim(min(beta_hat_1), max(beta_hat_1)) +
    labs(subtitle = paste0("sample size = ", N), x = NULL, y = "beta_hat_1")
}

point_plot <- plot_grid(plotlist = point_plot.list)
box_plot <- plot_grid(plotlist = box_plot.list)
print(point_plot)
print(box_plot)

##### sim2 #####

