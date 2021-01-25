rm(list = ls())
if (!require("MASS")) install.packages("MASS")
if (!require("dplyr")) install.packages("dplyr")
if (!require("plm")) install.packages("plm")

set.seed(123)

#######################
#####     DGP     #####
#######################

DGP3 <- function(N,T_,r){
  Lambda_i <- matrix(rnorm(n = r*N, mean = 0, sd = 1), nrow=r, ncol=N)
  F_t <- matrix(rnorm(n = r*T_, mean = 0, sd = 1), nrow=r, ncol=T_)
  Eps_it <- matrix(rnorm(n = T_*N, mean = 0, sd = 1), nrow=N, ncol=T_)
  X_it <- t(Lambda_i)%*%F_t+sqrt(r)*Eps_it
  X_list <- list()
  for (i in 1:N){
    X_list[[i]]<-X_it[i,]
  }
  return(X_list)
}

# Estimate F_tilde by X_list and r
F_tilde <- function(X_list, r){
  N <- length(X_list)
  T_ <- length(X_list[[1]])[1]
  XXT <- matrix(0, nrow=T_, ncol=T_)
  for (i in 1:N){
    X_i <- X_list[[i]]
    XXT <- XXT+ X_i %*% t(X_i)
  }
  eig <- eigen(XXT)
  F_tilde <- sqrt(T_)*eig$vectors[,order(eig$values, decreasing = T)[1:r]]
  return(F_tilde)
}

# Estimate lambda_tilde by F_tilde, X_list and r
lambda_tilde <- function(F_tilde, X_list, r){
  N <- length(X_list)
  T_ <- length(X_list[[1]])[1]
  Lambda_tilde <- matrix(0, nrow = N, ncol = r)
  for(i in 1:N){
    X_i <- X_list[[i]]
    Lambda_tilde[i,] <- t(F_tilde) %*% (X_i) / T_
  }
  return(Lambda_tilde)
}

# Define id for different g_functions
g <- function(N,T_,id){
  if(id==1){
    return((N+T_)/(N*T_)*log((N*T_)/(N+T_)))
  }
  if(id==2){
    return((N+T_)/(N*T_)*log(min(N,T_)))
  }
  if(id==3){
    return((log(min(N,T_)))/min(N,T_))
  }
  if(id==4){
    return(2/T_)
  }
  
  
}

# Define the sum of squared residuals function
VkF <- function(r,X_list,Lambda,F_){
  N <- length(X_list)
  T_ <- length(X_list[[1]])[1]
  
  X_it <- matrix(0, nrow=N, ncol=T_)
  for (i in (1:N)){
    X_it[i,]<-X_list[[i]]
  }
  return( sum((X_it-Lambda%*%t(F_))^2)/(N*T_))
}

# Define IC criteria 
IC <- function(r, X_list,id){
  N <- length(X_list)
  T_ <- length(X_list[[1]])[1]
  
  F_ <- F_tilde(X_list, r)
  Lambda <- lambda_tilde(F_, X_list, r)
  X_it <- matrix(0, nrow=N, ncol=T_)
  for (i in (1:N)){
    X_it[i,]<-X_list[[i]]
  }
  vkf <- VkF(r,X_list,Lambda,F_)
  IC_p <- log(vkf) +r*g(N,T_,id)
  return(IC_p)
  
}

# Define PC criteria 
PC <- function(r,X_list,rmax,id){
  N <- length(X_list)
  T_ <- length(X_list[[1]])[1]
  
  F_ <- F_tilde(X_list, r)
  Lambda <- lambda_tilde(F_, X_list, r)
  X_it <- matrix(0, nrow=N, ncol=T_)
  for (i in (1:N)){
    X_it[i,]<-X_list[[i]]
  }
  
  F__ <- F_tilde(X_list,rmax)
  Lambda_ <- lambda_tilde(F__,X_list,rmax)
  vkf <- VkF(r,X_list,Lambda,F_)
  sigma_sq <- VkF(rmax,X_list,Lambda_,F__)
  pc<- vkf+r*sigma_sq*g(N,T_,id)
  return(pc)
}

# Calculate r_hat by different criterias
r_hat <- function(rmax, X_list, panelty,id){
  v <- c()
  if(panelty=="IC"){
    for (r in (1:rmax)){
      v <- c(v, IC(r,X_list,id))
    }
    return(which.min(v))
  }
  if(panelty=="PC"){
    for (r in (1:rmax)){
      v <- c(v, PC(r,X_list,rmax,id))
    }
    return(which.min(v))
  }
}
  
# nsim<-1000
# r <- 2
# 
# for (i in 1:nsim){
#   r_hat_vector <- c()
#   N <- 100
#   T_ <- 60
#   X_list <- DGP3(N, T_,r)
#   rmax<-8
#   r_hat_vector <- c(r_hat_vector,r_hat(rmax, X_list,"PC",3))
# }
# 
# mean(r_hat_vector)

# Generate data frame to compare criterias
r <- 3
rmax <- 8
nsim <- 50
all_N <- c(100,100,200,500)
all_T <- c(40,60,60,60)
df <- data.frame(N=all_N,T_=all_T, PC1=NA,PC2=NA,PC3=NA,IC1=NA,IC2=NA,IC3=NA)
for(case in 1:length(all_N)){
  N <- all_N[case]
  T_ <- all_T[case]
  df$N[case] <- N
  df$T_[case] <- T_
  df_sim <- data.frame(PC1=rep(NA, nsim),PC2=NA,PC3=NA,IC1=NA,IC2=NA,IC3=NA)
  for(i in 1:nsim){
    X_list <- DGP3(N, T_,r)
    df_sim$PC1[i] <- r_hat(rmax, X_list,"PC",1)
    df_sim$PC2[i] <- r_hat(rmax, X_list,"PC",2)
    df_sim$PC3[i] <- r_hat(rmax, X_list,"PC",3)
    df_sim$IC1[i] <- r_hat(rmax, X_list,"IC",1)
    df_sim$IC2[i] <- r_hat(rmax, X_list,"IC",2)
    df_sim$IC3[i] <- r_hat(rmax, X_list,"IC",3)
  }
  df[case, 3:8] <- colMeans(df_sim)
}

df


## caculate r_hat

source("DGPs.R")
source("Methods.R")
source("Statistics.R")

beta_true <- c(1,3,5,2,4)
tolerance <- 0.0001
r0 <- 8
rmax <- 10
nsim<-50
T_ <- 100
N <- 100
model <- "model5"
r_hat_PC1_vector <- c()
r_hat_IC1_vector <- c()
for (i in 1:nsim) {
  sim_data <- DGP2(T_=T_, N=N, beta_true=beta_true, model)
  result_ls <- least_squares(sim_data$X_list, sim_data$Y_list, sim_data$df, tolerance, r0)
  beta_hat <- result_ls$beta_hat
  U_list <- mapply(function(Y_i, X_i) {Y_i - X_i %*% beta_hat}, sim_data$Y_list, sim_data$X_list, SIMPLIFY=F)
  
  r_hat_PC1_vector <- c(r_hat_PC1_vector, r_hat(rmax, U_list,"PC",1))
  r_hat_IC1_vector <- c(r_hat_IC1_vector, r_hat(rmax, U_list,"IC",1))
}
mean(r_hat_PC1_vector)
mean(r_hat_IC1_vector)
