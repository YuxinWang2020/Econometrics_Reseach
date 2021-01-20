rm(list = ls())

if (!require("MASS")) install.packages("MASS")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("akima")) install.packages("akima")
if (!require("plotly")) install.packages("plotly")
if (!require("reshape2")) install.packages("reshape2")
if (!require("plm")) install.packages("plm")


set.seed(123)
source("DGPs.R")
source("Methods.R")
source("Statistic.R")


##### Generate table #####
all_N <- c(100,100,100,100,10,20,50)
all_T <- c(10,20,50,100,100,100,100)
nsims <- 1000
beta_true <- c(1,3,5,2,4)
tolerance <- 0.0001
models <- c("model1","model2","model3","model4","model5")
sim_data_list <- as.list(rep(NA, length(models)))
names(sim_data_list) <- models
stat_ls_list <- stat_fe_list <- stat_ls_list <- table_ls_list <- table_fe_list <- sim_data_list
for(model in models){
  sim_data <- sim_dgp2_ls_fe(beta_true, tolerance, "model5", all_N, all_T, nsims)
  sim_data_list[[model]] <- sim_data
  stat_ls <- statistics(sim_data$df_beta_hat_ls, sim_data$df_sde, beta_true, all_N, all_T, nsims)
  stat_ls_list[[model]] <- stat_ls
  stat_fe <- statistics(sim_data$df_beta_hat_fe, NULL, beta_true, all_N, all_T, nsims, is.fe = T)
  stat_fe_list[[model]] <- stat_fe
  table_ls <- stat_ls$df_statistic[c("N","T_","mean.1","sde.1","mean.2","sde.2",
                                     "mean.3","sde.3","mean.4","sde.4","mean.5","sde.5")]
  colnames(table_ls) <- c("N","T","Mean β1=1","SD β1","Mean β2=3","SD β2","Mean μ=5","SD μ",
                          "Mean γ=2","SD γ","Mean =4","SD δ")
  table_ls_list[[model]] <- table_ls
  table_fe <- stat_fe$df_statistic[c("N","T_","mean.1","sde.1","mean.2","sde.2")]
  colnames(table_fe) <- c("N","T","Mean β1=1","SD β1","Mean β2=3","SD β2")
  table_fe_list[[model]] <- table_fe
}

