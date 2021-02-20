########################################
# Correlated Fixed Effects for Model 1 #
########################################

# model1 of DGP2 in the file "DGPs":
# y_it = beta1*x_it_1 + beta2*x_it_2 + alpha_i + eps     # p=2, X=(x1,x2)
# since model 1 doesn't contain iid varibles in fixed effects items, the result is same as the original model #

rm(list = ls())

if (!require("MASS")) install.packages("MASS")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("akima")) install.packages("akima")
if (!require("plotly")) install.packages("plotly")
if (!require("reshape2")) install.packages("reshape2")
if (!require("plm")) install.packages("plm")

# use JIT
library(compiler)
enableJIT(3)
setCompilerOptions(optimize=3)

set.seed(123)

### Data generating process ###
DGP <- function(T_, N, beta_true){
  p <- 2
  # Set parameters
  mu <- 0
  gamma <- 0
  delta <- 0
  iota <- rep(1, 2)
  mu1 <- mu2 <- c1 <- c2 <- 1
  
  # Generate variables
  Factor <- matrix(c(1,0), nrow=2, ncol=T_)
  Lambda <- rbind(matrix(rnorm(n = N, mean = 0, sd = 1), nrow=1, ncol=N), 1)
  
  Eta_1 <- matrix(rnorm(n = T_*N, mean = 0, sd = 1), nrow=N, ncol=T_)
  Eta_2 <- matrix(rnorm(n = T_*N, mean = 0, sd = 1), nrow=N, ncol=T_)
  Eps <- matrix(rnorm(n = T_*N, mean = 0, sd = 2), nrow=N, ncol=T_)
  e <- rnorm(n = N, mean = 0, sd = 1)
  eta <- rnorm(n = T_, mean = 0, sd = 1)
  # Calculate intermediate variables
  iota_Lambda <- crossprod(Lambda, iota)
  iota_Factor <- crossprod(iota, Factor)
  Lambda_Factor <- crossprod(Lambda, Factor)
  # Simulate data
  X_1 <- mu1 + c1 * Lambda_Factor + iota_Lambda %*% rep(1, T_) + rep(1, N) %*% iota_Factor + Eta_1
  X_2 <- mu2 + c2 * Lambda_Factor + iota_Lambda %*% rep(1, T_) + rep(1, N) %*% iota_Factor + Eta_2
  Y <- beta_true[1]*X_1 + beta_true[2]*X_2 + Lambda_Factor + Eps
  # Save all results to data frame
  X_df <- data.frame(x_it_1 = as.vector(t(X_1)),
                     x_it_2 = as.vector(t(X_2)))
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   y_it = as.vector(t(Y)),
                   X_df)
  # Save results to lists
  X_list <- list()
  Y_list <- list()
  for(i in 1:N){
    X_list[[i]] <- as.matrix(df[((i-1)*T_+1):(i*T_), 4:(3+p)])
    Y_list[[i]] <- as.matrix(df$y_it[((i-1)*T_+1):(i*T_)])
  }
  return(list(df=df, X_list=X_list, Y_list=Y_list))
}



### Interactive Fixed Effects Model ###

# same as the content in the file "Methods"
#Step 1:define funtion to calculate F_hat, dim of F_hat is (T_, r)
calculate_F_hat <- function(X_list, Y_list, beta_hat_0, r){
  N <- length(X_list)
  T_ <- dim(X_list[[1]])[1]
  
  WWT <- matrix(0, nrow=T_, ncol=T_)
  for(i in 1:N){
    X_i <- X_list[[i]]
    Y_i <- Y_list[[i]]
    W_i <- (Y_i - X_i %*% beta_hat_0)
    WWT <- WWT + W_i %*% t(W_i)
  }
  eig <- eigen(WWT)
  F_hat <- as.matrix(sqrt(T_)*eig$vectors[,order(eig$values, decreasing = T)[1:r]])
  return(F_hat)
}

#Step 2:define funtion to calculate Lambda_hat, dim of Lambda_hat is (r, N)
calculate_Lambda_hat <- function(X_list, Y_list, beta_hat_0, F_hat, r){
  N <- length(X_list)
  T_ <- dim(X_list[[1]])[1]
  p <- dim(X_list[[1]])[2]
  
  Lambda_hat <- matrix(NA, nrow = N, ncol = r)
  for(i in 1:N){
    X_i <- X_list[[i]]
    Y_i <- Y_list[[i]]
    Lambda_hat[i,] <- t(F_hat) %*% (Y_i - X_i %*% beta_hat_0) / T_
  }
  return(Lambda_hat)
}

#Step 3:define funtion to calculate Beta_hat
calculate_beta_hat <- function(X_list, Y_list, F_hat, Lambda_hat){
  N <- length(X_list)
  T_ <- dim(X_list[[1]])[1]
  p <- dim(X_list[[1]])[2]
  
  A <- matrix(0, nrow=p, ncol=p)
  B <- matrix(0, nrow=p, ncol=1)
  for(i in 1:N){
    X_i <- X_list[[i]]
    Y_i <- Y_list[[i]]
    A <- A + t(X_i) %*% X_i
    B <- B + t(X_i) %*% (Y_i - F_hat %*% Lambda_hat[i,])
  }
  beta_hat_0 <- solve(A) %*% B
  return(beta_hat_0)
}

#Step 4:calculate Beta_hat by iterations
least_squares <- function(X_list, Y_list, df, tolerance, r,beta_hat_0){
  # Initialize
  # p <- dim(X_list[[1]])[2]
  # formulate <- paste0("y_it ~ ", paste0("x_it_",c(1:p)) %>% paste(collapse = " + "), ifelse(p<=2, " + 0", ""))
  # beta_hat_0 <- plm(formulate, data=df, model="pooling")$coefficients %>% as.matrix()
  # if(p >=3){
  #   beta_hat_0[c(3,1,2)] <- beta_hat_0[c(1,2,3)]
  # }
  #beta_hat_list <- list(beta_hat_0)
  
  #beta_hat_0<-plm(y_it ~ x_it_1+x_it_2+x_it_0, data=df, model="pooling")$coefficients %>% as.matrix()
  
  e <- Inf
  while (e > tolerance) {
    F_hat <- calculate_F_hat(X_list, Y_list, beta_hat_0, r)
    Lambda_hat <- calculate_Lambda_hat(X_list, Y_list, beta_hat_0, F_hat, r)
    beta_hat <- calculate_beta_hat(X_list, Y_list, F_hat, Lambda_hat)
    #beta_hat_list[[length(beta_hat_list)+1]] <- beta_hat
    e <- norm(beta_hat - beta_hat_0, type = "F")
    beta_hat_0 <- beta_hat
  }
  
  return(list(beta_hat=beta_hat, F_hat=F_hat, Lambda_hat=Lambda_hat))
}

# calculate the mean value of beta_hat #
mean_value <- function(beta_hat_list){
  #Input: list of estimations
  #Output: vector of mean estimation, denoted by m
  N <- length(beta_hat_list) 
  p <- length(beta_hat_list[[1]])
  m <- rep(0, p)
  for (j in 1:p){
    for (i in 1:N){
      m[j] <- m[j]+beta_hat_list[[i]][j]
    }
  }
  return (1/N*m)
}

# calculate the rmse of beta_hat #
rmse <- function(beta_hat_list, real_beta){
  N <- length(beta_hat_list) 
  p <- length(real_beta)
  m <- rep(0,p)
  
  for (j in 1:p){
    for (i in 1:N){
      err_sqr<-(beta_hat_list[[i]][j]-real_beta[j])^2
      m[j]<- m[j] + err_sqr
    }
  }
  rmse<-m/N
  rmse<-sqrt(rmse)
  return(rmse)
}



### Generate result ###

r<-2 # Number of factors
nsim<-1000 # Number of simulations
tol<-0.001 # Iteration precision
real_beta<-c(1,3) # Regression coefficients

#Store the result in this data frame and set the names of the header
rst<-data.frame()

ncase<-7 # 7 combinations of different N & T
N_vec <- c(100,100,100,100,10,20,50) # Different sample sizes of N
T_vec <- c(10,20,50,100,100,100,100) # Different sample sizes of T

for(c in 1:ncase){
  N<-N_vec[c]
  T_<-T_vec[c]
  
  #list of interactive fix effects estimator
  beta_hat_list<-list()
  
  #list of pooled estimator
  beta_tilde_list<-list()
  
  for (i in 1:nsim){
    dgp<-DGP(T_,N, real_beta)
    X_list<-dgp$X_list
    Y_list<-dgp$Y_list
    df<-dgp$df
    
    
    beta_hat_0<-plm(y_it ~ x_it_1+x_it_2+0, data=df, model="pooling")$coefficients %>% as.matrix()
    
    #beta_hat_0<-plm(y_it ~ x_it_1+x_it_2+x_i+w_t+0, data=df, model="within")$coefficients %>% as.matrix()
    
    
    
    ls<-least_squares(X_list,Y_list,df,tol,r,beta_hat_0)
    beta_hat_list[[i]]<-ls$beta_hat
    beta_tilde_list[[i]]<-beta_hat_0
  }
  
  mv<-mean_value(beta_hat_list)
  sd<-rmse(beta_hat_list, real_beta)
  rowToAdd<-c(N,T_)
  for (i in 1:length(mv)){
    rowToAdd<-c(rowToAdd,mv[i])
    rowToAdd<-c(rowToAdd,sd[i])
  }
  rst<-rbind(rst,rowToAdd)
}

# Rename the data frame #
names(rst)<-c("N","T","beta1 mean", "beta1 sd", "beta2 mean", "beta2 sd")
rst

write.csv(rst, file = paste0( "model 1_", nsim, ".csv"), row.names = F)
