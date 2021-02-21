#####################################
#####     Factor Estimation     #####
#####################################

rm(list = ls())
if (!require("MASS")) install.packages("MASS")
if (!require("dplyr")) install.packages("dplyr")
if (!require("plm")) install.packages("plm")

set.seed(123)

# create dir
dir.create("../out", showWarnings = F)
dir.create("../out/tables", showWarnings = F)

source("DGPs.R")
source("Methods.R")


### Step 1: data generating process ###
# Y = beta*X + Lambda_i*Factor_t + eps
# U = Y - beta*X = Lambda_i*Factor_t + eps
DGP3 <- function(N,T_,r){
  Lambda_i <- matrix(rnorm(n = r*N, mean = 0, sd = 1), nrow=r, ncol=N)
  F_t <- matrix(rnorm(n = r*T_, mean = 0, sd = 1), nrow=r, ncol=T_)
  Eps_it <- matrix(rnorm(n = T_*N, mean = 0, sd = 1), nrow=N, ncol=T_)
  U_it <- t(Lambda_i)%*%F_t+sqrt(r)*Eps_it
  U_list <- list()
  for (i in 1:N){
    U_list[[i]]<-U_it[i,]
  }
  return(U_list)
}



### Step 2: calculation process following the paper Bai,Ng (2002) ###

# Estimate F_tilde by U_list and r #
F_tilde <- function(U_list, r){
  N <- length(U_list)
  T_ <- length(U_list[[1]])[1]
  XXT <- matrix(0, nrow=T_, ncol=T_)
  for (i in 1:N){
    U_i <- U_list[[i]]
    XXT <- XXT+ U_i %*% t(U_i)
  }
  eig <- eigen(XXT)
  F_tilde <- sqrt(T_)*eig$vectors[,order(eig$values, decreasing = T)[1:r]]
  return(F_tilde)
}

# Estimate lambda_tilde by F_tilde, U_list and r #
lambda_tilde <- function(F_tilde, U_list, r){
  N <- length(U_list)
  T_ <- length(U_list[[1]])[1]
  Lambda_tilde <- matrix(0, nrow = N, ncol = r)
  for(i in 1:N){
    U_i <- U_list[[i]]
    Lambda_tilde[i,] <- t(F_tilde) %*% (U_i) / T_
  }
  return(Lambda_tilde)
}

# Define id for different g_functions #
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

# Define the sum of squared residuals function #
VkF <- function(r,U_list,Lambda,F_){
  N <- length(U_list)
  T_ <- length(U_list[[1]])[1]
  
  U_it <- matrix(0, nrow=N, ncol=T_)
  for (i in (1:N)){
    U_it[i,]<-U_list[[i]]
  }
  return( sum((U_it-Lambda%*%t(F_))^2)/(N*T_))
}

# Define IC criteria #
IC <- function(r, U_list,id){
  N <- length(U_list)
  T_ <- length(U_list[[1]])[1]
  
  F_ <- F_tilde(U_list, r)
  Lambda <- lambda_tilde(F_, U_list, r)
  U_it <- matrix(0, nrow=N, ncol=T_)
  for (i in (1:N)){
    U_it[i,]<-U_list[[i]]
  }
  vkf <- VkF(r,U_list,Lambda,F_)
  IC_p <- log(vkf) +r*g(N,T_,id)
  return(IC_p)
  
}

# Define PC criteria #
PC <- function(r,U_list,rmax,id){
  N <- length(U_list)
  T_ <- length(U_list[[1]])[1]
  
  F_ <- F_tilde(U_list, r)
  Lambda <- lambda_tilde(F_, U_list, r)
  U_it <- matrix(0, nrow=N, ncol=T_)
  for (i in (1:N)){
    U_it[i,]<-U_list[[i]]
  }
  
  F__ <- F_tilde(U_list,rmax)
  Lambda_ <- lambda_tilde(F__,U_list,rmax)
  vkf <- VkF(r,U_list,Lambda,F_)
  sigma_sq <- VkF(rmax,U_list,Lambda_,F__)
  pc<- vkf+r*sigma_sq*g(N,T_,id)
  return(pc)
}

# Calculate r_hat by different criterias #
r_hat <- function(rmax, U_list, panelty,id){
  v <- c()
  if(panelty=="IC"){
    for (r in (1:rmax)){
      v <- c(v, IC(r,U_list,id))
    }
    return(which.min(v))
  }
  if(panelty=="PC"){
    for (r in (1:rmax)){
      v <- c(v, PC(r,U_list,rmax,id))
    }
    return(which.min(v))
  }
}



### Step 3: generate data frame to compare different criterias ###

# a). by using DGP3, we do a replication of Table 2 in Bai,Ng (2002), page 205 #

CompareCriterias <- function(r=3, rmax=8, nsim=1000, all_N=c(100,100,200,500,1000), all_T=c(40,60,60,60,60)){
  df <- data.frame(N=all_N,T_=all_T, PC1=NA,PC2=NA,PC3=NA,IC1=NA,IC2=NA,IC3=NA)
  for(case in 1:length(all_N)){
    N <- all_N[case]
    T_ <- all_T[case]
    df$N[case] <- N
    df$T_[case] <- T_
    df_sim <- data.frame(PC1=rep(NA, nsim),PC2=NA,PC3=NA,IC1=NA,IC2=NA,IC3=NA)
    for(i in 1:nsim){
      U_list <- DGP3(N, T_,r)
      df_sim$PC1[i] <- r_hat(rmax, U_list,"PC",1)
      df_sim$PC2[i] <- r_hat(rmax, U_list,"PC",2)
      df_sim$PC3[i] <- r_hat(rmax, U_list,"PC",3)
      df_sim$IC1[i] <- r_hat(rmax, U_list,"IC",1)
      df_sim$IC2[i] <- r_hat(rmax, U_list,"IC",2)
      df_sim$IC3[i] <- r_hat(rmax, U_list,"IC",3)
    }
    df[case, 3:8] <- colMeans(df_sim)
  }
  return(df)
}
df <- CompareCriterias()
df
write.csv(df, file = "../out/tables/determine_num_of_factors.csv", row.names = FALSE)

# b). by using DGP2, we use the same method for model4 in paper Bai(2009), page 1260 #
# y_it = mu + beta1*x_it_1 + beta2*x_it_2 + x_i*gamma + w_t*delta + Lambda_i*Factor_t + eps     
# p=5, X=(x1,x2,1,x_i,w_t)

CompareCriterias_DGP2 <- function(rmax=8,nsim=1000){
  # Set parameters #
  all_N <- c(100,100,100,100,10,20,50) # Different Sample sizes of N
  all_T <- c(10,20,50,100,100,100,100) # Different Sample sizes of T
  model="model4" # we use model 4 here
  beta_true=c(1,3,5,4,2)  # Regression coefficients
  tolerance=0.0001 # Iteration precision
  r0=8 # starting value of factor number
  
  # Initialize and loop #
  df <- data.frame(N=all_N,T_=all_T, PC1=NA,PC2=NA,PC3=NA,IC1=NA,IC2=NA,IC3=NA)
  for(case in 1:length(all_N)){
    N <- all_N[case]
    T_ <- all_T[case]
    df$N[case] <- N
    df$T_[case] <- T_
    df_sim <- data.frame(PC1=rep(NA, nsim),PC2=NA,PC3=NA,IC1=NA,IC2=NA,IC3=NA)
    for(i in 1:nsim){
      sim_data <- DGP2(T_=T_, N=N, beta_true=beta_true, model)
      result_ls <- least_squares(sim_data$X_list, sim_data$Y_list, sim_data$df, tolerance, r0, model)
      beta_hat <- result_ls$beta_hat
      U_list <- mapply(function(Y_i, X_i) {Y_i - X_i %*% beta_hat}, sim_data$Y_list, sim_data$X_list, SIMPLIFY=F)
      
      df_sim$PC1[i] <- r_hat(rmax, U_list,"PC",1)
      df_sim$PC2[i] <- r_hat(rmax, U_list,"PC",2)
      df_sim$PC3[i] <- r_hat(rmax, U_list,"PC",3)
      df_sim$IC1[i] <- r_hat(rmax, U_list,"IC",1)
      df_sim$IC2[i] <- r_hat(rmax, U_list,"IC",2)
      df_sim$IC3[i] <- r_hat(rmax, U_list,"IC",3)
    }
    df[case, 3:8] <- colMeans(df_sim)
  }
  return(df)
}
df_dgp2 <- CompareCriterias_DGP2()
df_dgp2
write.csv(df_dgp2, file = "../out/tables/determine_num_of_factors_dgp2.csv", row.names = FALSE)



### Step 4: save results to a file ###
save.image(file = "../out/tables/factorEstimation.RData")
