##########################
#####     Statistics     #####
##########################

#####  Compute mean squared error  #####
#input: data frame or matrix of estimations and real parameters
#output: mean squared error
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

rmse <- function(est_df, real_para){
  mse <- rep(0, 4)
  N <- nrow(est_df)
  real_para <- real_para[1:ncol(est_df)]
  residual <- t(t(est_df) - real_para) # est_df substract real_para by row
  rmse <- colMeans(residual^2)
  return(rmse)
}

#####  Compute standard error  #####
#Page 1240, define funtion to calculate a_ik
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

#Page 1241, define funtion to calculate Zi
calculate_Z <- function(X_list, N, M, a){
  Z_list <- list()
  for(i in 1:N){
    # dims of Z_i is c(T, p)
    Z_list[[i]] <- M %*% X_list[[i]] - 1/N * M %*% X_list[[i]] * sum(a[i,])
  }
  return(Z_list)
}

#Page 1246 & 1252, define funtion to calculate D0 and D1
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

#Page 1246, define funtion to calculate covariance matrix of beta_hat
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

#####  Simulation  #####

sim_dgp1_fe <- function(beta_true, all_N, all_T, nsims){
  # Data frame to save every beta_hat
  p <- length(beta_true)
  df_beta_hat_fe <- data.frame(T_ = rep(all_T, each = nsims),
                               N = rep(all_N, each = nsims),
                               sim = rep(1:nsims, times = length(all_N)),
                               beta = matrix(NA,ncol=p))
  # Loop over all_N and all_T and c(1:nsims) for simulation
  loop_count <- 1
  for(case in 1:length(all_N)){
    N <- all_N[case]
    T_ <- all_T[case]
    for(h in 1:nsims){
      sim_data <- DGP1(T_=T_, N=N, beta_true=beta_true)
      result_fe <- OLS_FE(sim_data$X_list, sim_data$Y_list)
      df_beta_hat_fe[loop_count, 4:(3+p)] <- result_fe$beta_hat
      loop_count <- loop_count + 1
    }
  }
  return(list(df_beta_hat_fe=df_beta_hat_fe))
}

sim_dgp2_ls_fe <- function(beta_true, tolerance, r, model, all_N, all_T, nsims, need.sde, need.fe){
  p <- ifelse(model == "model4", 5, ifelse(model == "model3", 3, 2))
  # Data frame to save every beta_hat
  df_beta_hat_ls <- data.frame(T_ = rep(all_T, each = nsims),
                            N = rep(all_N, each = nsims),
                            sim = rep(1:nsims, times = length(all_N)),
                            beta = matrix(NA,ncol=p))
  df_beta_hat_fe <- data.frame(T_ = rep(all_T, each = nsims),
                            N = rep(all_N, each = nsims),
                            sim = rep(1:nsims, times = length(all_N)),
                            beta = matrix(NA,ncol=p))
  df_sde <- data.frame(T_ = rep(all_T, each = nsims),
                       N = rep(all_N, each = nsims),
                       sim = rep(1:nsims, times = length(all_N)),
                       sde = matrix(NA,ncol=p))
  # Loop over all_N and all_T and c(1:nsims) for simulation
  loop_count <- 1
  for(case in 1:length(all_N)){
    N <- all_N[case]
    T_ <- all_T[case]
    
    for(h in 1:nsims){
      sim_data <- DGP2(T_=T_, N=N, beta_true=beta_true, model)
      result_ls <- least_squares(sim_data$X_list, sim_data$Y_list, sim_data$df, tolerance, r)
      df_beta_hat_ls[loop_count, 4:(3+p)] <- result_ls$beta_hat
      if(need.fe){
        # X_list_no_singular <- lapply(sim_data$X_list, function(X_i) X_i[,1:2])
        # df_beta_hat_fe[loop_count, 4:5] <- OLS_FE(X_list_no_singular, sim_data$Y_list)$beta_hat
        df_no_singular <- sim_data$df[,1:5]
        df_beta_hat_fe[loop_count, 4:5] <- OLS_FE2(df_no_singular)$beta_hat
      }
      if(need.sde){
        ls_sde <- calculate_sde(sim_data$X_list, sim_data$Y_list, result_ls$beta_hat, result_ls$F_hat, result_ls$Lambda_hat)
        df_sde[loop_count, 4:(3+p)] <- sqrt(diag(ls_sde))
      }
      loop_count <- loop_count + 1
    }
  }
  return(list(df_beta_hat_ls=df_beta_hat_ls, df_sde=df_sde, df_beta_hat_fe=df_beta_hat_fe))
}

#####  Statistics  #####
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
