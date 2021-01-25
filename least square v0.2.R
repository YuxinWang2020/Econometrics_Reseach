rm(list = ls())


if (!require("MASS")) install.packages("MASS")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("akima")) install.packages("akima")
if (!require("plotly")) install.packages("plotly")
if (!require("reshape2")) install.packages("reshape2")
if (!require("plm")) install.packages("plm")
library(RSpectra)
library(psych)
set.seed(123)

################################
#####     Data Process     #####
################################

#Turn a list to a 3-dimensional matrix
listToMatrix <- function(X_list){
  N <- length(X_list)
  k <- length(X_list[[1]][1,])
  T_ <- length(X_list[[1]])/k
  A <- array(0, c(T_, N, k))
  for(i in 1:N){
    A[,i,]<-X_list[[i]]
  }
  return(A)
}

#######################
#####     DGP     #####
#######################


DGP <- function(T_, N, beta_true){
  # Set parameters
  r<-2  #number of factors
  K<-length(beta_true) #number of parameters
  mu = 0 #grand mean
  gamma <- 0
  delta <- 0
  iota <- rep(1, r)
  mu1 <- mu2 <- c1 <- c2 <- 1
  # Generate variables
  F_t <- matrix(rnorm(n = r*T_, mean = 0, sd = 1), nrow=r, ncol=T_)
  Lambda_i <- matrix(rnorm(n = r*N, mean = 0, sd = 1), nrow=r, ncol=N)
  Eta_it_1 <- matrix(rnorm(n = T_*N, mean = 0, sd = 1), nrow=N, ncol=T_)
  Eta_it_2 <- matrix(rnorm(n = T_*N, mean = 0, sd = 1), nrow=N, ncol=T_)
  Eps_it <- matrix(rnorm(n = T_*N, mean = 0, sd = 4), nrow=N, ncol=T_)
  e_i <- rnorm(n = N, mean = 0, sd = 1)
  eta_t <- rnorm(n = T_, mean = 0, sd = 1)
  # Calculate intermediate variables
  iota_Lambda_i = crossprod(iota, Lambda_i)
  iota_F_t = crossprod(iota, F_t)
  dim(iota_Lambda_i) <- c(N, 1)
  dim(iota_F_t) <- c(1, T_)
  iota_Lambda_it <- iota_Lambda_i %*% rep(1, T_)
  iota_F_it <- rep(1, N) %*% iota_F_t
  x_i <- iota_Lambda_i + e_i
  w_t <- iota_F_t + eta_t
  gamma_x_it <- gamma * x_i %*% rep(1, T_)
  delta_w_it <- delta * rep(1, N) %*% w_t
  Lambda_F_it <- crossprod(Lambda_i, F_t)
  # Simulate data
  X_it_1 <- mu1 + c1 * Lambda_F_it + iota_Lambda_it + iota_F_it + Eta_it_1
  X_it_2 <- mu2 + c2 * Lambda_F_it + iota_Lambda_it + iota_F_it + Eta_it_2
  Y_it <- beta_true[1]*X_it_1 + beta_true[2]*X_it_2 + mu + gamma_x_it + delta_w_it + 
    Lambda_F_it + Eps_it
  # Save all results to data frame
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   y_it = as.vector(t(Y_it)),
                   #x_it_0 = 1,
                   x_it_1 = as.vector(t(X_it_1)),
                   x_it_2 = as.vector(t(X_it_2)))
  
  # Save results to lists
  X_list <- list()
  Y_list <- list()
  
  for(i in 1:N){
    X_list[[i]] <- as.matrix(df[((i-1)*T_+1):(i*T_), 4:(3+K)])
    Y_list[[i]] <- as.matrix(df$y_it[((i-1)*T_+1):(i*T_)])
  }
  
  X <- listToMatrix(X_list)
  Y <- matrix(unlist(Y_list), nrow=T_, ncol=N)
  
  return(list(df=df, X=X, Y=Y, X_list=X_list, Y_list=Y_list, F_t=F_t, Lambda_i=Lambda_i))
}


#########################
####   Estimation    ####
#########################

#Step 1: Compute XXT

xxinv <- function(X){
  T_ <- dim(X)[1]
  N <- dim(X)[2]
  p <- dim(X)[3]
  xx <- matrix(0, p, p)
  if (p==1){
    xx[1,1] <- tr(t(X)%*%X)
  } else {
    for (k in 1:p){
      X1 <- X[,,k]
      for (m in k:p){
        X2=X[,,m]
        xx[k,m] <- tr(t(X1)%*%X2)
        if (k<m){
          xx[m,k]<-xx[k,m]
        }
      }
    }
  }
  return(solve(xx))
}

#Step 2: Calculate F and Lambda
#Calculate F and lambda
calculate_F_lambda <- function(X, r){
  T_ <- dim(X)[1]
  N <- dim(X)[2]
  
  if (T_<N){
    XX=X%*%t(X)/(N*T_)
    SVD <- svd(XX)
    Factor <- SVD$u[,1:r]*sqrt(T_)
    VNT <- SVD$d
    lambda <- t(X)%*%Factor/T_
  } else{
    XX=t(X)%*%X/(N*T_)
    SVD <- svd(XX)
    lambda <- SVD$u[,1:r]*sqrt(N)
    VNT <- SVD$d
    Factor <- X%*%lambda/N
  }
  return(list(f=Factor, l=lambda, ev=VNT))
}

#Step 3: Calculate beta
calculate_beta_hat <- function(X, XXinv, y, f, l){
  T_ <- dim(X)[1]
  N <- dim(X)[2]
  p <- dim(X)[3]
  
  xy<-matrix(0, ncol=1, nrow=p)
  
  for(k in 1:p){
    xy[k] <- tr(t(X[,,k])%*%(y-f%*%t(l)))
  }
  return(XXinv%*%xy)
}


#Step 4: Iteration of step 1 to 3
least_squares <- function(X, Y, tolerance, beta_hat_old, r){
  
  beta_hat_list <- list(beta_hat_old)
  A <- xxinv(X)
  e <- Inf
  while (e > tolerance) {
    U <- Y
    for(k in 1:2){
      U <- U-X[,,k]*beta_hat_old[k]
    }
    
    FL <- calculate_F_lambda(U, r)
    
    beta_hat <- calculate_beta_hat(X, A, Y, FL$f, FL$l)
    
    beta_hat_list[[length(beta_hat_list)+1]] <- beta_hat
    e <- norm(beta_hat - beta_hat_old, type = "F")
    beta_hat_old <- beta_hat
  }
  return(list(beta_hat=beta_hat, beta_hat_list=beta_hat_list))
}


###################################
#####     Result Analysis     #####
###################################

#####  Compute mean squared error  #####


mse <- function(est_list, real_para){
  #input: list of estimations and real parameters
  #output: mean squared error vector mse. mse[1] means the mse of the first parameter  
  #N is set to be the number of parameters
  N <- length(est_list) 
  
  #Initialize 
  mse <- rep(0,length(real_para))
  
  #For each estimation i, calculate diff, the difference between the estimation and the true value
  #For each parameter j, calculate the mse and store it in mse[j]
  for (i in 1:N){
    
    diff <- est_list[[i]]-real_para
    for (j in 1:length(real_para)){
      mse[j] <- mse[j]+diff[j]^2
    }
  }
  mse <- 1/N*mse
  return(mse)
}


mean_value <- function(est_list){
  #Input: list of estimations
  #Output: vector of mean estimation, denoted by m
  N <- length(est_list) 
  k <- length(est_list[[1]])
  m <- rep(0, k)
  for (j in 1:k){
    for (i in 1:N){
      m[j] <- m[j]+est_list[[i]][j]
    }
  }
  return (1/N*m)
}



#############################################
#####     Start of the main program     #####
#############################################

#Interactive fixed effect

T_ <- 100
N <- 100
tol <-0.0001


#List that stores the estimate of beta_hat in each regression
beta_hat_list <- list()

#number of regressions
n_reg <- 1


for (i in 1:n_reg){
  dgp<-DGP(T_, N, c(1,3))
  X <- dgp$X
  Y <- dgp$Y
  df<-dgp$df
  
  
  
  beta_hat_old <- plm(y_it ~ x_it_1+x_it_2+0, data=df, model="pooling")$coefficients %>% as.matrix()
  beta_hat_list[[i]]<-least_squares(X,Y,0.0001,beta_hat_old, 2)$beta_hat
  
  
  
}

mean_value(beta_hat_list)