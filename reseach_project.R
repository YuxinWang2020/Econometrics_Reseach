# K=5 #
rm(list = ls())

if (!require("MASS")) install.packages("MASS")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("akima")) install.packages("akima")
if (!require("plotly")) install.packages("plotly")
if (!require("reshape2")) install.packages("reshape2")
if (!require("plm")) install.packages("plm")
# if (!require("RSpectra")) install.packages("RSpectra")
# if (!require("psych")) install.packages("psych")
# library(RSpectra)
# library(psych)

set.seed(123)

#######################
#####     DGP     #####
#######################

DGP1 <- function(T_, N, beta_true){
  # Set parameters
  p <- 3
  range_x <- c(0,1)
  alpha_i <- 1
  mu_u_i <- rep(0, T_)
  sigma_u_i <- diag(rep(1, T_))
  # Simulate data 
  X <- lapply(as.list(1:N), function(i) array(runif(p*T_,range_x[1],range_x[2]), dim = c(T_, p)))
  U <- lapply(as.list(1:N), function(i) mvrnorm(n = 1, mu = mu_u_i, Sigma = sigma_u_i))
  Y <- mapply(function(X_i, u_i) alpha_i + X_i %*% beta_true + u_i, X, U, SIMPLIFY=F)
  # Alternatively, we can write in the following way:
  # X_i <- array(runif(p*T_,range_x[1],range_x[2]), dim = c(T_, p))
  # u_i <- mvrnorm(n = 1, mu = mu_u_i, Sigma = sigma_u_i)
  # Y_i <- alpha_i + X_i %*% beta_true + u_i
  # Save all results to data frame
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   y_it = NA) %>% cbind(matrix(NA,nrow=N*T_,ncol=p))
  colnames(df)[4:(3+p)] <- paste0("x_it_",c(0:(p-1)))
  # Loop over N & T_
  for(i in 1:N){
    X_i <- X[[i]]
    u_i <- U[[i]]
    Y_i <- Y[[i]]
    for(t in 1:T_){
      df[(i-1)*T_ + t, 4:(3+p)] <- X_i[t,]
      df$y_it[(i-1)*T_ + t] <- Y_i[t]
    }
  }
  return(list(X_list = X, Y_list = Y, df=df))
}

DGP2 <- function(T_, N, beta_true){
  if(T_ < 2 || N < 2)
    stop("T or N less than 2")
  alpha_i <- runif(N, min = 0, max = 1)
  eps_t <- rnorm(T_, mean = 0, sd = 1)
  z0 <- 0
  z_t <- c(0.5*z0 + eps_t[1])
  for(i in 2:T_)
    z_t[i] <- 0.5*z_t[i-1] + eps_t[i]
  v_i <- rnorm(N, mean = 0, sd = 2)
  e_t <- rnorm(T_ + 1, mean = 0, sd = 1)
  mu <- 5
  w_t <- mu + e_t[2:(T_+1)] + e_t[1:T_] / 3
  x0 <- 1
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   y_it = NA,
                   x_it_0 = x0,
                   x_it_1 = NA)
  for(i in 1:N){
    for(t in 1:T_){
      x_it <- sqrt(alpha_i[i]) * z_t[t]
      df$x_it_1[(i-1)*T_ + t] <- x_it
      u_it <- v_i[i] + w_t[t]
      df$y_it[(i-1)*T_ + t] <- alpha_i[i] + beta_true[1] * x0 + beta_true[2] * x_it + u_it
    }
  }
  # Save results to lists
  X_list <- list()
  Y_list <- list()
  p <- 2
  for(i in 1:N){
    X_list[[i]] <- as.matrix(df[((i-1)*T_+1):(i*T_), 4:(3+p)])
    Y_list[[i]] <- as.matrix(df$y_it[((i-1)*T_+1):(i*T_)])
  }
  return(list(df=df, X_list=X_list, Y_list=Y_list))
}

DGP3 <- function(T_, N, beta_true=c(1,3,5,2,4), p=5){
  # Set parameters
  r<-2
  mu <- beta_true[3]
  gamma <- beta_true[4]
  delta <- beta_true[5]
  iota <- rep(1, r)
  mu1 <- mu2 <- c1 <- c2 <- 1
  # Generate variables
  F_t <- matrix(rnorm(n = r*T_, mean = 0, sd = 1), nrow=r, ncol=T_)
  Lambda_i <- matrix(rnorm(n = r*N, mean = 0, sd = 1), nrow=r, ncol=N)
  Eta_it_1 <- matrix(rnorm(n = T_*N, mean = 0, sd = 1), nrow=N, ncol=T_)
  Eta_it_2 <- matrix(rnorm(n = T_*N, mean = 0, sd = 1), nrow=N, ncol=T_)
  Eps_it <- matrix(rnorm(n = T_*N, mean = 0, sd = 2), nrow=N, ncol=T_)
  e_i <- rnorm(n = N, mean = 0, sd = 1)
  eta_t <- rnorm(n = T_, mean = 0, sd = 1)
  # Calculate intermediate variables
  iota_Lambda_i <- crossprod(iota, Lambda_i)
  iota_F_t <- crossprod(iota, F_t)
  dim(iota_Lambda_i) <- c(N, 1)
  dim(iota_F_t) <- c(1, T_)
  iota_Lambda_it <- iota_Lambda_i %*% rep(1, T_)
  iota_F_it <- rep(1, N) %*% iota_F_t
  x_i <- iota_Lambda_i + e_i
  w_t <- iota_F_t + eta_t
  Lambda_F_it <- crossprod(Lambda_i, F_t)
  # Simulate data
  X_it_1 <- mu1 + c1 * Lambda_F_it + iota_Lambda_it + iota_F_it + Eta_it_1
  X_it_2 <- mu2 + c2 * Lambda_F_it + iota_Lambda_it + iota_F_it + Eta_it_2
  X_it_3 <- matrix(1, nrow=N, ncol=T_)
  X_it_4 <- x_i %*% rep(1, T_)
  X_it_5 <- rep(1, N) %*% w_t
  Y_it <- beta_true[1]*X_it_1 + beta_true[2]*X_it_2 + mu*X_it_3 + gamma*X_it_4 + delta*X_it_5 + 
    Lambda_F_it + Eps_it
  # Save all results to data frame
  X_df <- data.frame(x_it_1 = as.vector(t(X_it_1)),
                     x_it_2 = as.vector(t(X_it_2)),
                     x_it_3 = as.vector(t(X_it_3)),
                     x_it_4 = as.vector(t(X_it_4)),
                     x_it_5 = as.vector(t(X_it_5)))
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   y_it = as.vector(t(Y_it)),
                   X_df[,1:p])
  # Save results to lists
  X_list <- list()
  Y_list <- list()
  for(i in 1:N){
    X_list[[i]] <- as.matrix(df[((i-1)*T_+1):(i*T_), 4:(3+p)])
    Y_list[[i]] <- as.matrix(df$y_it[((i-1)*T_+1):(i*T_)])
  }
  return(list(df=df, X_list=X_list, Y_list=Y_list))
}

#####################
####   Method    ####
#####################

#####  Introduction   #####

### Fixed Effects Model ###
OLS_FE <- function(df){
  # Initialize
  N <- length(unique(df$i))
  T_ <- length(unique(df$t))
  p <- ncol(df) - 3
  A <- array(0, dim = c(p,p))
  B <- rep(0, p)
  # Loop over N & T_
  for(i in 1:N){
    X_i <- as.matrix(df[((i-1)*T_+1):(i*T_), 4:(3+p)])
    Y_i <- as.matrix(df$y_it[((i-1)*T_+1):(i*T_)])
    x_i_mean <- colMeans(X_i)
    y_i_mean <- mean(Y_i)
    for(t in 1:T_){
      A <- A + (X_i[t,] - x_i_mean) %*% t(X_i[t,] - x_i_mean)
      B <- B + (X_i[t,] - x_i_mean) * (Y_i[t] - y_i_mean)
    }
  }
  beta_hat_fe <- solve(A) %*% B
  return(list(beta_hat = beta_hat_fe))
}

# alternatively, we can use plm package for FE estimation
OLS_FE2 <- function(df){
  p <- ncol(df) - 3
  data <- pdata.frame(df,index=c("i","t"))
  formulate <- reformulate(response = c("y_it"), termlabels =
                             paste0("x_it_",c(0:(p-1))) %>% paste(collapse = " + "))
  result <- plm(formulate, data=data,
                effect = "individual",model="within")
  return(list(beta_hat = as.matrix(result$coefficients)))
}


#####  Interacctive Fixed Effect Methods  #####

### Least Squares Model ###
#Step 1:define funtion to calculate F_hat, dim of F_hat is (T_, r)
calculate_F_hat <- function(X_list, Y_list, beta_hat, r){
  N <- length(X_list)
  T_ <- dim(X_list[[1]])[1]
  p <- dim(X_list[[1]])[2]
  
  WWT <- matrix(0, nrow=T_, ncol=T_)
  for(i in 1:N){
    X_i <- X_list[[i]]
    Y_i <- Y_list[[i]]
    W_i <- (Y_i - X_i %*% beta_hat)
    WWT <- WWT + W_i %*% t(W_i)
  }
  eig <- eigen(WWT)
  F_hat <- sqrt(T_)*eig$vectors[,order(eig$values, decreasing = T)[1:r]]
  return(F_hat)
}

#Step 2:define funtion to calculate Lambda_hat, dim of Lambda_hat is (r, N)
calculate_Lambda_hat <- function(X_list, Y_list, beta_hat, F_hat, r){
  N <- length(X_list)
  T_ <- dim(X_list[[1]])[1]
  p <- dim(X_list[[1]])[2]
  
  Lambda_hat <- matrix(NA, nrow = N, ncol = r)
  for(i in 1:N){
    X_i <- X_list[[i]]
    Y_i <- Y_list[[i]]
    Lambda_hat[i,] <- t(F_hat) %*% (Y_i - X_i %*% beta_hat) / T_
  }
  return(Lambda_hat)
}

#Step 3:define funtion to calculate Beta_hat
calculate_beta_hat <- function(X_list, Y_list, F_, Lambda){
  N <- length(X_list)
  T_ <- dim(X_list[[1]])[1]
  p <- dim(X_list[[1]])[2]
  
  A <- matrix(0, nrow=p, ncol=p)
  B <- matrix(0, nrow=p, ncol=1)
  for(i in 1:N){
    X_i <- X_list[[i]]
    Y_i <- Y_list[[i]]
    A <- A + t(X_i) %*% X_i
    B <- B + t(X_i) %*% (Y_i - F_ %*% Lambda[i,])
  }
  beta_hat <- solve(A) %*% B
  return(beta_hat)
}

#Step 5:calculate Beta_hat by iterations
least_squares <- function(X_list, Y_list, df, tolerance, r){
  # Initialize
  p <- dim(X_list[[1]])[2]
  formulate <- paste0("y_it ~ ", paste0("x_it_",c(1:p)) %>% paste(collapse = " + "), ifelse(p<=2, " + 0", ""))
  beta_hat_0 <- plm(formulate, data=df, model="pooling")$coefficients %>% as.matrix()
  if(p >=3){
    beta_hat_0[c(3,1,2)] <- beta_hat_0[c(1,2,3)]
  }
  beta_hat_list <- list(beta_hat_0)
  e <- Inf
  while (e > tolerance) {
    F_hat <- calculate_F_hat(X_list, Y_list, beta_hat_0, r)
    Lambda_hat <- calculate_Lambda_hat(X_list, Y_list, beta_hat_0, F_hat, r)
    beta_hat <- calculate_beta_hat(X_list, Y_list, F_hat, Lambda_hat)
    beta_hat_list[[length(beta_hat_list)+1]] <- beta_hat
    e <- norm(beta_hat - beta_hat_0, type = "F")
    beta_hat_0 <- beta_hat
  }
  
  return(list(beta_hat=beta_hat, beta_hat_list=beta_hat_list, F_hat=F_hat, Lambda_hat=Lambda_hat))
}


##########################
#####     Result     #####
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
sim_dgp3_ls <- function(beta_true, p, tolerance, r, all_N, all_T, nsims){
  # Data frame to save every beta_hat
  df_beta_hat <- data.frame(T_ = rep(all_T, each = nsims),
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
      sim_data <- DGP3(T_=T_, N=N, beta_true=beta_true)
      result <- least_squares(sim_data$X_list, sim_data$Y_list, sim_data$df, tolerance, r)
      beta_hat <- result$beta_hat
      df_beta_hat[loop_count, 4:(3+p)] <- beta_hat
      
      ls_sde <- calculate_sde(sim_data$X_list, sim_data$Y_list, beta_hat, result$F_hat, result$Lambda_hat)
      df_sde[loop_count, 4:(3+p)] <- sqrt(diag(ls_sde))
      loop_count <- loop_count + 1
    }
  }
  return(list(df_beta_hat=df_beta_hat, df_sde=df_sde))
}

#####  Statistics  #####
statistics <- function(df_beta_hat, df_sde, beta_true, all_N, all_T, nsims){
  # Initialize
  p <- ncol(df_beta_hat) - 3
  # Data frame to save statistics variable
  df_statistic <- data.frame(T_ = all_T,
                             N = all_N,
                             mean = matrix(ncol = p),
                             bias = matrix(ncol = p),
                             sde = matrix(ncol = p),
                             ci_l = matrix(ncol = p),
                             ci_u = matrix(ncol = p),
                             mse = NA)
  findcol <- function(column){
    return(substr(colnames(df_statistic),1,nchar(column))==column)
  }
  loop_count <- 1
  for(case in 1:length(all_N)){
    N <- all_N[case]
    T_ <- all_T[case]
    row_range <- c((loop_count*nsims-nsims+1) : (loop_count*nsims)) # rows range of beta_hat for N and T_
    
    df_statistic[loop_count, findcol("mean")] <- colMeans(df_beta_hat[row_range,4:(3+p)], na.rm = T)
    df_statistic[loop_count, findcol("bias")] <- 
      abs(df_statistic[loop_count, findcol("mean")] - beta_true[1:p]) / beta_true[1:p]
    df_statistic$mse[loop_count] <- mse(df_beta_hat[row_range, 4:(3+p)], beta_true)
    df_statistic[loop_count, findcol("sde")] <- colMeans(df_sde[row_range,4:(3+p)], na.rm = T)
    df_statistic[loop_count, findcol("ci_l")] <- df_statistic[loop_count, findcol("mean")] -
      qnorm(0.975,0,1)*df_statistic[loop_count, findcol("sde")]
    df_statistic[loop_count, findcol("ci_u")] <- df_statistic[loop_count, findcol("mean")] +
      qnorm(0.975,0,1)*df_statistic[loop_count, findcol("sde")]
    loop_count <- loop_count + 1
  }
  return(list(df_statistic = df_statistic))
}

#####  Test  #####
dgp2 <- DGP2(10, 100, beta_true = c(1,2))
OLS_FE2(dgp2$df)

dgp1 <- DGP1(T_ = 3,
             N = 400,  beta_true = c(1,2,3))
OLS_FE(dgp1$df)$beta_hat
OLS_FE2(dgp1$df)$beta_hat
all.equal(OLS_FE(dgp1$df)$beta_hat, OLS_FE2(dgp1$df)$beta_hat)

T_ <- 30
N <- 50
tol <-0.0001
beta_true <- c(1,3,5,2,4)
p <- 5
df_beta_hat <- data.frame(matrix(nrow=0, ncol = 5)) #List that stores the estimate of beta_hat in each regression
n_reg <- 10 #number of regressions
for (i in 1:n_reg){
  dgp<-DGP3(T_, N, beta_true, p)
  df<-dgp$df
  X_list <- dgp$X_list
  Y_list <- dgp$Y_list
  ls <- least_squares(X_list, Y_list, df, tol, r=10)
  df_beta_hat[i, 1:p] <- ls$beta_hat
  
  F_hat <- ls$F_hat
  Lambda_hat <- ls$Lambda_hat
  ls_sde <- calculate_sde(X_list, Y_list, ls$beta_hat, F_hat, Lambda_hat)
}
(MSE <- mse(df_beta_hat, beta_true))
(MEAN <- colMeans(df_beta_hat))


##### Generate table #####
all_N <- c(100,100,100,100,10,20,50)
all_T <- c(10,20,50,100,100,100,100)
nsims <- 5
beta_true <- c(1,3,5,2,4)
p <- 5
r <- 4
tolerance <- 0.0001
sim1 <- sim_dgp3_ls(beta_true, p, tolerance, r, all_N, all_T, nsims)
stat1 <- statistics(sim1$df_beta_hat, sim1$df_sde, beta_true, all_N, all_T, nsims)
View(stat1$df_statistic)
table1 <- stat1$df_statistic[c("N","T_","mean.1","sde.1","mean.2","sde.2",
                               "mean.3","sde.3","mean.4","sde.4","mean.5","sde.5")]
colnames(table1) <- c("N","T","Mean β1=1","SD β1","Mean β2=3","SD β2","Mean μ=5","SD μ",
                      "Mean γ=2","SD γ","Mean =4","SD δ")



###########################
####   Visualization   ####
###########################

##### sim1 #####
all_n <- c(1000, 2000, 3000)
nsims <- 100
beta_true <- c(1:3)
sim1 <- sim(p = 3, beta_true = beta_true, all_n = all_n, all_T = c(3), nsims = nsims)
beta_hat_1 <- sim1$beta_hat_1
beta_true_1 <- beta_true[1]
beta_hat_1_df <- data.frame(t(beta_hat_1[,1,]))
colnames(beta_hat_1_df) <- all_n
beta_hat_1_df <- melt(beta_hat_1_df, measure.vars=colnames(beta_hat_1_df), variable.name = "N")
# generate point and box plot for every N
plot.heigth <- max(beta_true_1 - min(beta_hat_1_df$value), max(beta_hat_1_df$value) - beta_true_1)
point_plot <- ggplot(data = beta_hat_1_df) +
  geom_point(aes(y = value, x = c(1:(nsims*length(all_n))), color = N, shape = N)) +
  geom_hline(yintercept = beta_true_1, color = I("black")) +
  ylim(beta_true_1 - plot.heigth, beta_true_1 + plot.heigth) +
  labs(title = "point plot for every N", x = "iteration", y = "beta_hat_1") +
  theme(axis.text.x = element_blank())
print(point_plot)

box_plot <- ggplot(data = beta_hat_1_df) +
  geom_boxplot(aes(y = value, x = N, color = N)) +
  ylim(min(beta_hat_1_df$value), max(beta_hat_1_df$value)) +
  labs(title = "box plot for every N", y = "beta_hat_1")
print(box_plot)

##### sim2 #####
sim2 <- sim(p = 3, beta_true = c(1:3), all_n = seq(1000, 10000, 1000), all_T = seq(3,30,3), nsims = 1)
beta_hat_1 <- sim2$beta_hat_1
beta_true_1 <- 1
summary_beta_hat_1 <- sim2$summary_beta_hat_1

im <- with(summary_beta_hat_1, interp(T_,N,bias))
with(im,image(x,y,z, xlab = "T_", ylab = "N"))

plot_ly(x=im$x, y=im$y, z=im$z, type="surface") %>% add_surface() %>% layout(scene = list(
  xaxis = list(title = 'T'),
  yaxis = list(title = 'N'),
  zaxis = list(title = 'bias')))

