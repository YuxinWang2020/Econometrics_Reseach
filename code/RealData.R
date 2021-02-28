#################################
##    Real Data Application    ##
#################################

rm(list = ls())

if (!require("plm")) install.packages("plm")
library(dplyr)
source("Methods.R")

### Data Management ###

data("Cigar", package = "plm")

# Data cleaning #
Cigar_clean <- pdata.frame(Cigar, index = c("state", "year"))

Cigar_clean$log_C_it <- log(Cigar_clean$sales * Cigar_clean$pop / Cigar_clean$pop16)
Cigar_clean$log_C_it_lag1 <- plm::lag(Cigar_clean$log_C_it)  # Log C_{i,t-1}
Cigar_clean$log_P_it <- log(Cigar_clean$price / Cigar_clean$cpi)
Cigar_clean$log_Y_it <- log(Cigar_clean$ndi / Cigar_clean$cpi)
Cigar_clean$log_Pn_it <- log(Cigar_clean$pimin / Cigar_clean$cpi)

Cigar_clean <- na.omit(Cigar_clean) # remove NA


### Pooled Estimation ###
# use plm package for estimation

ols_estimate <- plm("log_C_it ~ log_C_it_lag1 + log_P_it + log_Pn_it + log_Y_it",
    Cigar_clean, model="pooling")
ols_estimate$coefficients

within_estimate <- plm("log_C_it ~ log_C_it_lag1 + log_P_it + log_Pn_it + log_Y_it",
    Cigar_clean, model="within",effect="twoways")
within_estimate$coefficients


### Interactive-effects Estimation ###
# use model 3 of DGP 2 in the file "DGPs" for estimation
least_squares <- function(X_list, Y_list, df, tolerance, r, beta_hat_0){
  # Initialize
  beta_hat_list <- list(beta_hat_0)
  e <- Inf
  
  while (e > tolerance) {
    F_hat <- calculate_F_hat(X_list, Y_list, beta_hat_0, r)
    Lambda_hat <- calculate_Lambda_hat(X_list, Y_list, beta_hat_0, F_hat, r)
    beta_hat <- calculate_beta_hat(X_list, Y_list, F_hat, Lambda_hat)
    beta_hat_list[[length(beta_hat_list)+1]] <- beta_hat
    e <- norm(beta_hat - beta_hat_0, type = "2")
    beta_hat_0 <- beta_hat
  }
  
  return(list(beta_hat=beta_hat, beta_hat_list=beta_hat_list, F_hat=F_hat, Lambda_hat=Lambda_hat))
}

# Set parameters #
tolerance = 0.001
states <- levels(index(Cigar_clean)$state) # unique state in Cigar_clean

## with log_C_it_lag1 ##
df <- mutate(Cigar_clean, Intercept = 1) # add constant	term for interactive effects estimate
X_list <- list()
Y_list <- list()
# Loop states. Store rows with same state in one matrix and put it in the list.
for(i in 1:length(states)){
  X_list[[i]] <- as.matrix( df[df$state==states[i], c("Intercept","log_C_it_lag1","log_P_it","log_Pn_it","log_Y_it")] )
  Y_list[[i]] <- as.matrix( df$log_C_it[df$state==states[i]] )
}
beta_hat_0 <- plm("log_C_it ~ log_C_it_lag1 + log_P_it + log_Pn_it + log_Y_it",
                  data=df, model="pooling")$coefficients %>% as.matrix()
interactive_effects_estimate <- least_squares(X_list, Y_list, df, tolerance, r = 2, beta_hat_0)
interactive_effects_estimate$beta_hat


### Factor Estimation ###

# Set parameters #
tolerance = 0.001
states <- levels(index(Cigar_clean)$state) # unique state in Cigar_clean

## with log_C_it_lag1 ##
df <- mutate(Cigar_clean, Intercept = 1) # add constant	term for interactive effects estimate
X_list <- list()
Y_list <- list()

# Loop states. Store rows with same state in one matrix and put it in the list.
for(i in 1:length(states)){
  X_list[[i]] <- as.matrix( df[df$state==states[i], c("Intercept","log_C_it_lag1","log_P_it","log_Pn_it","log_Y_it")] )
  Y_list[[i]] <- as.matrix( df$log_C_it[df$state==states[i]] )
}
beta_hat_0 <- plm("log_C_it ~ log_C_it_lag1 + log_P_it + log_Pn_it + log_Y_it",
                  data=df, model="pooling")$coefficients %>% as.matrix()

rs <- c(1:10)
df_factor_estimation <- data.frame(r = c(0, rs), Intercept=NA, log_C_it_lag1=NA, log_P_it=NA, log_Pn_it=NA, log_Y_it=NA,
                                   row.names = c(0, rs))
df_factor_estimation[1,2:6] <- beta_hat_0
for(r in rs){
  interactive_effects_estimate <- least_squares(X_list, Y_list, df, tolerance, r = r, beta_hat_0)
  df_factor_estimation[r+1,2:6] <- interactive_effects_estimate$beta_hat
}
df_factor_estimation



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



### Compare Different Criterias ###

# Set parameters #
r0s <- c(1:10) # starting value of factor number
rmax <- 12

# Initialize and loop #
df_r_hat <- data.frame(r0=r0s, PC1=NA,PC2=NA,PC3=NA,IC1=NA,IC2=NA,IC3=NA)
for(i in 1:length(r0s)){
  beta_hat <- t(as.matrix(df_factor_estimation[df_factor_estimation$r==r0s[i], 2:6]))
  U_list <- mapply(function(Y_i, X_i) {Y_i - X_i %*% beta_hat}, Y_list, X_list, SIMPLIFY=F)
  
  df_r_hat$PC1[i] <- r_hat(rmax, U_list,"PC",1)
  df_r_hat$PC2[i] <- r_hat(rmax, U_list,"PC",2)
  df_r_hat$PC3[i] <- r_hat(rmax, U_list,"PC",3)
  df_r_hat$IC1[i] <- r_hat(rmax, U_list,"IC",1)
  df_r_hat$IC2[i] <- r_hat(rmax, U_list,"IC",2)
  df_r_hat$IC3[i] <- r_hat(rmax, U_list,"IC",3)
}

df_r_hat
