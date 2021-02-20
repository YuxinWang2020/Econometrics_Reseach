####################################
# Change Starting Value of Model 2 #
####################################

# model2: y_it = beta1*x_it_1 + beta2*x_it_2 + alpha_i + z_t + eps     
# p=2, X=(x1,x2)

rm(list = ls())

if (!require("dplyr")) install.packages("dplyr")
if (!require("MASS")) install.packages("MASS")
if (!require("plm")) install.packages("plm")

# set.seed(123)
source("DGPs.R")
source("Methods.R")
source("Statistics.R")

### Interactive Fixed Effects Model ###

### ls1: Starting value is pooling estimator ###

least_squares1 <- function(X_list, Y_list, df, tolerance, r, model){
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
    e <- norm(beta_hat - beta_hat_0, type = "2")
    beta_hat_0 <- beta_hat
  }
  
  return(list(beta_hat=beta_hat, beta_hat_list=beta_hat_list, F_hat=F_hat, Lambda_hat=Lambda_hat))
}

### ls2: Starting value is twoways estimator ### 
least_squares2 <- function(X_list, Y_list, df, tolerance, r, model){
  # Initialize
  p <- dim(X_list[[1]])[2]
  formulate <- paste0("y_it ~ ", paste0("x_it_",c(1:p)) %>% paste(collapse = " + "), ifelse(p<=2, " + 0", ""))
  beta_hat_0 <- plm(formulate, data=df,effect = "twoways", model="within")$coefficients %>% as.matrix()
  
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
    e <- norm(beta_hat - beta_hat_0, type = "2")
    beta_hat_0 <- beta_hat
  }
  
  return(list(beta_hat=beta_hat, beta_hat_list=beta_hat_list, F_hat=F_hat, Lambda_hat=Lambda_hat))
}

# Set parameters for small tests #
T_<-100 # Sample size of T
N<-100 # Sample size of N
r<-2 # Number of factors 
nsims<- 10 # Number of simulations
tolerance<-0.0001 # Iteration precision
model <- "model2" # we use model 2 in chaning starting values
beta_true<-c(1,3,5,0,0) # Regression coefficients

beta_hat_df <- data.frame(method = rep(as.factor(c("ls1","ls2","fe3")), nsims), beta_hat = matrix( NA, ncol = 3))
# "ls1","ls2": starting values from pooling and twoways estimator respectively, for interactive fixed effect model
# "fe3": starting values from twoways estimator, for simple fixed effects model

for (i in 1:nsims){
  dgp<-DGP2(T_,N, beta_true, model)
  X_list<-dgp$X_list
  Y_list<-dgp$Y_list
  df<-dgp$df
  
  ls1 <- least_squares1(X_list,Y_list,df,tolerance,r, model)
  ls2 <- least_squares2(X_list,Y_list,df,tolerance,r, model)
  fe3 <- OLS_FE3(df)
  beta_hat_df[i*3 - 2, 2:(length(ls1$beta_hat)+1)]<-ls1$beta_hat
  beta_hat_df[i*3 - 1, 2:(length(ls2$beta_hat)+1)]<-ls2$beta_hat
  beta_hat_df[i*3, 2:(length(fe3$beta_hat)+1)]<-fe3$beta_hat
}
View(beta_hat_df)
beta_hat_df %>% group_by(method) %>% summarise(mean(beta_hat.1), mean(beta_hat.2))


### 3: Generate table for comparison ###

# Set parameters #
all_N <- c(100,100,100,100,10,20,50) # Different sample sizes of N
all_T <- c(10,20,50,100,100,100,100) # Different sample sizes of T
nsims<-1000 # Number of simulations
r<-2 # Number of factors 
tolerance<-0.0001 # Iteration precision
model <- "model2" # we use model 2 in chaning starting values
beta_true<-c(1,3,5,0,0) # Regression coefficients

# Loop over all_N and all_T and c(1:nsims) for simulation
df_beta_hat_ls <- data.frame(method = rep(as.factor(c("ls1","ls2","fe3")), nsims*length(all_N)),
                             T_ = rep(all_T, each = nsims*3),
                             N = rep(all_N, each = nsims*3),
                             sim = rep(1:nsims, times=length(all_N), each=3),
                             beta = matrix(NA,ncol=3))
loop_count <- 1
for(case in 1:length(all_N)){
  N <- all_N[case]
  T_ <- all_T[case]
  for(h in 1:nsims){
    dgp<-DGP2(T_,N, beta_true, model)
    X_list<-dgp$X_list
    Y_list<-dgp$Y_list
    df<-dgp$df
    
    ls1 <- least_squares1(X_list,Y_list,df,tolerance,r, model)
    ls2 <- least_squares2(X_list,Y_list,df,tolerance,r, model)
    fe3 <- OLS_FE3(df)
    df_beta_hat_ls[loop_count*3 - 2, 5:(length(ls1$beta_hat)+4)]<-ls1$beta_hat
    df_beta_hat_ls[loop_count*3 - 1, 5:(length(ls2$beta_hat)+4)]<-ls2$beta_hat
    df_beta_hat_ls[loop_count*3, 5:(length(fe3$beta_hat)+4)]<-fe3$beta_hat
    loop_count <- loop_count + 1
  }
}
sv <- df_beta_hat_ls %>% group_by(N,T_,method) %>% summarise(mean(beta.1, na.rm=T), mean(beta.2, na.rm=T))
View(sv)

write.csv(sv, file = "test_starting_value/sv.csv", row.names = FALSE)
