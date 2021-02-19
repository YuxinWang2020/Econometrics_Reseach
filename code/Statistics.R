##############################
#####     Statistics     #####
##############################

library(foreach)
if (!require("doParallel")) install.packages("doParallel")

#####  Compute mean squared error  #####

# input: data frame or matrix of estimations and real parameters
# output: mean squared error
mse <- function(est_df, real_para){
  mse <- 0
  N <- nrow(est_df)
  real_para <- real_para[1:ncol(est_df)]
  for (i in 1:N){
    mse <- mse + norm(est_df[i,]-real_para, type="2")^2
  }
  mse <- 1/N * mse
  return(mse)
}

# calculate the mean of root-mean-square errors of different simulations
rmse <- function(est_df, real_para){
  mse <- rep(0, 4)
  N <- nrow(est_df)
  real_para <- real_para[1:ncol(est_df)]
  residual <- t(t(est_df) - real_para) # est_df substract real_para by row
  rmse <- colMeans(residual^2)
  rmse <- sqrt(rmse)
  return(rmse)
}



#####  Compute standard error of interactive estimator #####

#Page 1240, define funtion to calculate a_ik.
calculate_a <- function(N, Lambda_hat){
  A <- solve(crossprod(Lambda_hat, Lambda_hat) / N)
  a <- matrix(NA, nrow = N, ncol = N)
  for(i in 1:N){
    for(p in 1:N){
      a[i,p] <- t(Lambda_hat[i,]) %*% A %*% Lambda_hat[p,]
    }
  }
  return(a)
}

#Page 1241, define funtion to calculate Zi.
calculate_Z <- function(X_list, N, M, a){
  Z_list <- list()
  for(i in 1:N){
    # dims of Z_i is c(T, p)
    Z_list[[i]] <- M %*% X_list[[i]] - 1/N * M %*% X_list[[i]] * sum(a[i,])
  }
  return(Z_list)
}

#Page 1246 & 1252, define funtion to calculate D0 and D1.
calculate_D <- function(X_list, Y_list, N, T_, p, beta_hat, F_hat, Lambda_hat, Z){
  sita_square <- rep(0, N)
  for(i in 1:N){
    for(t in 1:T_){
      sita_square[i]  <- sita_square[i] + 1/T_ * (Y_list[[i]][t] - crossprod(beta_hat, X_list[[i]][t,])
                                                  - crossprod(Lambda_hat[i,], F_hat[t,])) ^ 2
    }
  }
  
  D0 <- matrix(0, nrow = p, ncol = p)
  D1 <- matrix(0, nrow = p, ncol = p)
  for(i in 1:N){
    for(t in 1:T_){
      A <- 1/(N*T_) * tcrossprod(Z[[i]][t,], Z[[i]][t,])
      D0 <- D0 + A
      D1 <- D1 + sita_square[i] * A
    }
  }
  return(list(D0=D0, D1=D1))
}

#Page 1246, define funtion to calculate covariance matrix of beta_hat.
calculate_sde <- function(X_list, Y_list, beta_hat, F_hat, Lambda_hat){
  N <- length(X_list)
  T_ <- dim(X_list[[1]])[1]
  p <- dim(X_list[[1]])[2]
  
  M <- diag(1, nrow = T_) - tcrossprod(F_hat, F_hat) / T_
  a <- calculate_a(N, Lambda_hat)
  Z <- calculate_Z(X_list, N, M, a)
  D <- calculate_D(X_list, Y_list, N, T_, p, beta_hat, F_hat, Lambda_hat, Z)
  sde <- solve(D$D0) %*% D$D1 %*% solve(t(D$D0))
  return(sde)
}

# calculate the mean value of beta_hat
mean_value <- function(est_list){
  #Input: list of estimations
  #Output: vector of mean estimation, denoted by m
  N <- length(est_list) 
  p <- length(est_list[[1]])
  m <- rep(0, p)
  for (j in 1:p){
    for (i in 1:N){
      m[j] <- m[j]+est_list[[i]][j]
    }
  }
  return (1/N*m)
}




#####  Monte Carlo Simulations  #####

# simulations for DGP1
sim_dgp1_fe <- function(beta_true, all_N, all_T, nsims){
  # init for parallel computing
  cl<-makeCluster(detectCores()) # create cluster with all cpu cores
  registerDoParallel(cl) # register cluster for DoParallel
  func <- function(case){ # function to run in cluster
    T_ <- T_N_sim$T_[case]
    N <- T_N_sim$N[case]
    sim_data <- DGP1(T_=T_, N=N, beta_true=beta_true)
    result_fe <- OLS_FE(sim_data$X_list, sim_data$Y_list)
    return(unlist(c(T_N_sim[case,], result_fe$beta_hat))) # combine results into a vector
  }
  clusterExport(cl=cl, varlist=c("DGP1", "OLS_FE"), envir=environment()) # export vars to cluster environment
  
  # Loop over all_N and all_T and c(1:nsims) for simulation
  T_N_sim <- data.frame(T_ = rep(all_T, each=nsims),
                        N = rep(all_N, each=nsims),
                        sim = rep(1:nsims, times=length(all_T)))
  # Loop over T_N_sim, use `rbind` to combine results from every simulation
  df_beta_hat_fe <- foreach(case=1:nrow(T_N_sim), .combine='rbind') %dopar%
    func(case) %>% as.data.frame() # columns is : T_, N, sim, beta.1, beta.2, beta...
  
  stopImplicitCluster() # stop cluster
  colnames(df_beta_hat_fe) <- c("T_", "N", "sim", paste0("beta.", 1:length(beta_true)))
  return(list(df_beta_hat_fe=df_beta_hat_fe))
}

# simulations for DGP2
sim_dgp2_ls_fe <- function(beta_true, tolerance, r, model, all_N, all_T, nsims, need.sde, need.fe){
  # init for parallel computing
  cl<-makeCluster(detectCores()) # create cluster with all cpu cores
  registerDoParallel(cl) # register cluster for DoParallel
  func <- function(case){ # function to run in cluster
    T_ <- T_N_sim$T_[case]
    N <- T_N_sim$N[case]
    sim_data <- DGP2(T_=T_, N=N, beta_true=beta_true, model)
    result_ls <- least_squares(sim_data$X_list, sim_data$Y_list, sim_data$df, tolerance, r, model)
    beta_hat_ls <- result_ls$beta_hat
    beta_hat_fe <- rep(NA, 2)
    if(need.fe){
      df_no_singular <- sim_data$df[,1:5] # select T,N,y,x1,x2 from sim_data$df
      # if model1, use OLS_FE2
      if(model == "model1"){
        fe_result <- OLS_FE2(df_no_singular)
        beta_hat_fe <- fe_result$beta_hat
      } else {
        fe_result <- OLS_FE2(df_no_singular)
        beta_hat_fe <- fe_result$beta_hat
      }
    }
    ls_sde <- rep(NA, p)
    if(need.sde){
      ls_sde <- calculate_sde(sim_data$X_list, sim_data$Y_list, result_ls$beta_hat, result_ls$F_hat, result_ls$Lambda_hat)
      ls_sde <- sqrt(diag(ls_sde))
    }
    return(unlist(c(T_N_sim[case,], beta_hat_ls, beta_hat_fe, ls_sde)))
  }
  clusterExport(cl=cl, varlist=c(ls(.GlobalEnv)), envir=environment()) # export all vars to cluster environment
  
  p <- ifelse(model == "model4", 5, ifelse(model == "model3", 3, 2)) # p is length of X
  # Loop over all_N and all_T and c(1:nsims) for simulation
  T_N_sim <- data.frame(T_ = rep(all_T, each=nsims),
                        N = rep(all_N, each=nsims),
                        sim = rep(1:nsims, times=length(all_T)))
  # Loop over T_N_sim, use `rbind` to combine results from every simulation
  beta.ls_beta.fe_ls.sde <- foreach(case=1:nrow(T_N_sim), .combine='rbind', .packages=c("plm", "dplyr")) %dopar%
    func(case) %>% as.data.frame() # columns is : T_, N, sim, beta_hat_ls(length=p), beta_hat_fe(length=2), ls_sde(length=p)
  stopImplicitCluster() # stop cluster
  
  df_beta_hat_ls <- beta.ls_beta.fe_ls.sde[, c(1:3, 4:(3+p))]
  df_beta_hat_fe <- beta.ls_beta.fe_ls.sde[, c(1:3, (4+p):(5+p))]
  df_sde <- beta.ls_beta.fe_ls.sde[, c(1:3, (6+p):(5+2*p))]
  colnames(df_beta_hat_ls) <- c("T_", "N", "sim", paste0("beta.", 1:p))
  colnames(df_beta_hat_fe) <- c("T_", "N", "sim", paste0("beta.", 1:2))
  colnames(df_sde) <- c("T_", "N", "sim", paste0("sde.", 1:p))
  return(list(df_beta_hat_ls=df_beta_hat_ls, df_sde=df_sde, df_beta_hat_fe=df_beta_hat_fe))
}




#####  Statistics  #####
# generate statistics of each N & T, take the mean of different simulations, and store them in a data frame
# including mean, bias, rmse, standard error and cofidence interval

statistics <- function(df_beta_hat, df_sde, beta_true, all_N, all_T, nsims){
  # Initialize
  p <- ncol(df_beta_hat) - 3
  beta_true <- beta_true[1:p]
  # Data frame to save statistics variable
  df_statistic <- data.frame(T_ = all_T,
                             N = all_N,
                             mean = matrix(ncol = p),
                             bias = matrix(ncol = p),
                             sde = matrix(ncol = p),
                             ci_l = matrix(ncol = p),
                             ci_u = matrix(ncol = p),
                             rmse = matrix(ncol = p))
  findcol <- function(column){
    return(substr(colnames(df_statistic),1,nchar(column))==column)
  }
  loop_count <- 1
  for(case in 1:length(all_N)){
    N <- all_N[case]
    T_ <- all_T[case]
    row_range <- c((loop_count*nsims-nsims+1) : (loop_count*nsims)) # rows range of beta_hat for N and T_
    
    # df_statistic[loop_count, paste0("mean.", 1:p)] <- colMeans(df_beta_hat[row_range,4:(3+p)], na.rm = T)
    
    df_statistic[loop_count, findcol("mean")] <- colMeans(df_beta_hat[row_range,4:(3+p)], na.rm = T)
    df_statistic[loop_count, findcol("bias")] <- 
      t(abs(t(df_statistic[loop_count, findcol("mean")]) - beta_true) / beta_true)
    if(!is.null(df_sde)){
      df_statistic[loop_count, findcol("sde")] <- colMeans(df_sde[row_range,4:(3+p)], na.rm = T)
      df_statistic[loop_count, findcol("ci_l")] <- df_statistic[loop_count, findcol("mean")] -
        qnorm(0.975,0,1)*df_statistic[loop_count, findcol("sde")]
      df_statistic[loop_count, findcol("ci_u")] <- df_statistic[loop_count, findcol("mean")] +
        qnorm(0.975,0,1)*df_statistic[loop_count, findcol("sde")]
    }
    df_statistic[loop_count, findcol("rmse")] <- rmse(df_beta_hat[row_range, 4:(3+p)], beta_true)
    loop_count <- loop_count + 1
  }
  return(df_statistic)
}
