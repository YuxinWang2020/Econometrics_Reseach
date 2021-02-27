###################################
#####     Generate Tables     #####
###################################

rm(list = ls())

if (!require("dplyr")) install.packages("dplyr")
if (!require("MASS")) install.packages("MASS")
if (!require("plm")) install.packages("plm")

# use JIT
library(compiler)
setCompilerOptions(optimize=3)
enableJIT(3)

# create dir
dir.create("../out", showWarnings = F)
dir.create("../out/tables", showWarnings = F)

set.seed(123)
source("DGPs.R")
source("Methods.R")
source("Statistics.R")


##### Generate table #####
coefficients <- c("beta1", "beta2", "mu", "gamma", "delta")


### 1: Repeatition of Table 1 in paper Bai(2009), page 1261 ###

# Set parameters #
all_N <- c(100,100,100,100,10,20,50) # Different Sample sizes of N
all_T <- c(10,20,50,100,100,100,100) # Different Sample sizes of T
nsims <- 1000 # Number of simulations
beta_true <- c(1,3,5,2,4) # Regression coefficients
tolerance <- 0.0001 # Iteration precision
r <- 2 # Number of factors
models <- c("model2") # Different models defined in the file DGPs

# Initialize #
sim_data_list1 <- as.list(rep(NA, length(models)))
names(sim_data_list1) <- models
table_loop_models_list <- stat_ls_list1 <- stat_fe_list1 <- sim_data_list1
select_statistics <- list(colName = c("mean", "rmse"), presentName = c("Mean", "SD"))

# Loop over different models & N & T, for both within estimator & interactive fe estimator #
for(model in models){
  sim_data <- sim_dgp2_ls_fe(beta_true, tolerance, r, model, all_N, all_T, nsims, need.sde=F, need.fe=T) # set need.sde=T if sde is contained in select_statistics
  sim_data_list1[[model]] <- sim_data

  stat_ls <- statistics(sim_data$df_beta_hat_ls, sim_data$df_sde, beta_true, all_N, all_T, nsims)
  stat_ls_list1[[model]] <- stat_ls
  stat_fe <- statistics(sim_data$df_beta_hat_fe, NULL, beta_true, all_N, all_T, nsims)
  stat_fe_list1[[model]] <- stat_fe

  p <- ifelse(model == "model4", 5, ifelse(model == "model3", 3, 2))
  table_ls <- stat_ls[c("N","T_", paste0(rep(select_statistics$colName, p), ".",
                                         rep(1:p, each=length(select_statistics$colName))))]
  colnames(table_ls) <- c("N","T_", paste(rep(select_statistics$presentName, p),
                                          rep(coefficients[1:p], each=length(select_statistics$presentName))))
  p <- 2
  table_fe <- stat_fe[c("N","T_", paste0(rep(select_statistics$colName, p), ".",
                                         rep(1:p, each=length(select_statistics$colName))))]
  colnames(table_fe) <- c("N","T_", paste(rep(select_statistics$presentName, p),
                                          rep(coefficients[1:p], each=length(select_statistics$presentName))))

  table_loop_models <- bind_rows(cbind(method="Least Squares", table_ls),
                             cbind(method="Fixed Effects", table_fe))
  write.csv(table_loop_models, file = paste0("../out/tables/table_",model,".csv"), row.names = FALSE)
  table_loop_models_list[[model]] <- table_loop_models
}

# 
# 
# ### 2: Change number of factors and compare the statistical results ###
# 
# # Set parameters #
# r_N <- c(50) # Sample sizes of N
# r_T <- c(50) # Sample sizes of T
# nsims <- 1000 # Number of simulations
# beta_true <- c(1,3,5,2,4) # Regression coefficients
# tolerance <- 0.0001 # Iteration precision
# model <- "model4" # we use model 4 in chaning r
# p <- ifelse(model == "model4", 5, ifelse(model == "model3", 3, 2)) #
# rs <- c(1:10) # range of r, from 1 to 10
# 
# # Initialize #
# sim_data_list2 <- as.list(rep(NA, length(rs)))
# names(sim_data_list2) <- paste0("r",rs)
# stat_ls_list2 <- sim_data_list2
# select_statistics <- list(colName = c("mean", "rmse", "sde"), presentName = c("Mean", "SD", "SDE"))
# table_loop_r <- data.frame(matrix(NA, nrow = 0, ncol = 0))
# 
# # Loop over different factor numbers #
# for(i in 1:length(rs)){
#   r <- rs[i]
#   sim_data <- sim_dgp2_ls_fe(beta_true, tolerance, r, model, r_N, r_T, nsims, need.sde=T, need.fe=F)
#   sim_data_list2[[i]] <- sim_data
#   stat_ls <- statistics(sim_data$df_beta_hat_ls, sim_data$df_sde, beta_true, r_N, r_T, nsims)
#   stat_ls_list2[[i]] <- stat_ls
#   
#   table_ls <- stat_ls[c("N","T_", paste0(rep(select_statistics$colName, p), ".",
#                                          rep(1:p, each=length(select_statistics$colName))))]
#   colnames(table_ls) <- c("N","T_", paste(rep(select_statistics$presentName, p),
#                                           rep(coefficients[1:p], each=length(select_statistics$presentName))))
#   table_loop_r <- bind_rows(table_loop_r,  cbind(r=r, table_ls))
# }
# write.csv(table_loop_r, file = "../out/tables/table_loop_r.csv", row.names = FALSE)
# 


### 3: Save results to a file###
save.image(file = "../out/tables/tables.RData")
