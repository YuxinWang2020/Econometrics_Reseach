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
source("Statistics.R")


##### Generate table #####

# sim1: Loop models & N & T for ls and fe
all_N <- c(100,100,100,100,10,20,50)
all_T <- c(10,20,50,100,100,100,100)
nsims <- 10
beta_true <- c(1,3,5,2,4)
tolerance <- 0.0001
r <- 2
models <- c("model1","model2","model3","model4","model5")
sim_data_list1 <- as.list(rep(NA, length(models)))
names(sim_data_list1) <- models
stat_ls_list1 <- stat_fe_list1 <- sim_data_list1
for(model in models){
  sim_data <- sim_dgp2_ls_fe(beta_true, tolerance, r, model, all_N, all_T, nsims, need_sde=T, need_fe=T)
  sim_data_list1[[model]] <- sim_data
  stat_ls <- statistics(sim_data$df_beta_hat_ls, sim_data$df_sde, beta_true, all_N, all_T, nsims, is.fe = F)
  stat_ls_list1[[model]] <- stat_ls
  stat_fe <- statistics(sim_data$df_beta_hat_fe, NULL, beta_true, all_N, all_T, nsims, is.fe = T)
  stat_fe_list1[[model]] <- stat_fe
}
table_ls_model1 <- stat_ls_list1[["model1"]][c("N","T_","mean.1","sde.1","mean.2","sde.2")]
colnames(table_ls_model1) <- c("N","T","Mean β1=1","SD β1","Mean β2=3","SD β2")

table_ls_model2 <- stat_ls_list1[["model2"]][c("N","T_","mean.1","sde.1","mean.2","sde.2", "mean.3","sde.3")]
colnames(table_ls_model2) <- c("N","T","Mean β1=1","SD β1","Mean β2=3","SD β2","Mean μ=5","SD μ")

table_ls_model3 <- stat_ls_list1[["model3"]][c("N","T_","mean.1","sde.1","mean.2","sde.2", "mean.3","sde.3")]
colnames(table_ls_model3) <- c("N","T","Mean β1=1","SD β1","Mean β2=3","SD β2","Mean μ=5","SD μ")

table_ls_model4 <- stat_ls_list1[["model4"]][c("N","T_","mean.1","sde.1","mean.2","sde.2", "mean.3","sde.3")]
colnames(table_ls_model4) <- c("N","T","Mean β1=1","SD β1","Mean β2=3","SD β2","Mean μ=5","SD μ")

table_ls_model5 <- stat_ls_list1[["model5"]][c("N","T_","mean.1","sde.1","mean.2","sde.2",
                                              "mean.3","sde.3","mean.4","sde.4","mean.5","sde.5")]
colnames(table_ls_model5) <- c("N","T","Mean β1=1","SD β1","Mean β2=3","SD β2","Mean μ=5","SD μ",
                               "Mean γ=2","SD γ","Mean =4","SD δ")

table_fe_model1 <- stat_fe_list1[["model1"]][c("N","T_","mean.1","sde.1","mean.2","sde.2")]
colnames(table_fe_model1) <- c("N","T","Mean β1=1","SD β1","Mean β2=3","SD β2")

table_fe_model2 <- stat_fe_list1[["model2"]][c("N","T_","mean.1","sde.1","mean.2","sde.2")]
colnames(table_fe_model2) <- c("N","T","Mean β1=1","SD β1","Mean β2=3","SD β2")

table_fe_model3 <- stat_fe_list1[["model3"]][c("N","T_","mean.1","sde.1","mean.2","sde.2")]
colnames(table_fe_model3) <- c("N","T","Mean β1=1","SD β1","Mean β2=3","SD β2")

table_fe_model4 <- stat_fe_list1[["model4"]][c("N","T_","mean.1","sde.1","mean.2","sde.2")]
colnames(table_fe_model4) <- c("N","T","Mean β1=1","SD β1","Mean β2=3","SD β2")

table_fe_model5 <- stat_fe_list1[["model5"]][c("N","T_","mean.1","sde.1","mean.2","sde.2")]
colnames(table_fe_model5) <- c("N","T","Mean β1=1","SD β1","Mean β2=3","SD β2")

################
## Loop for r ##
################
all_N <- c(50)
all_T <- c(50)
nsims <- 5
beta_true <- c(1,3,5,2,4)
tolerance <- 0.0001
model <- "model5"
rs <- c(2:10)
sim_data_list2 <- as.list(rep(NA, length(rs)))
names(sim_data_list2) <- paste0("r",rs)
table_loop_r_list <- stat_ls_list2 <- sim_data_list2
for(i in 1:length(rs)){
  r <- rs[i]
  sim_data <- sim_dgp2_ls_fe(beta_true, tolerance, r, model, all_N, all_T, nsims, need_sde=T, need_fe=F)
  sim_data_list2[[i]] <- sim_data
  stat_ls <- statistics(sim_data$df_beta_hat_ls, sim_data$df_sde, beta_true, all_N, all_T, nsims, is.fe = F)
  stat_ls_list2[[i]] <- stat_ls
  table_loop_r <- stat_ls[c("N","T_","mean.1","sde.1","mean.2","sde.2",
                                                 "mean.3","sde.3","mean.4","sde.4","mean.5","sde.5")]
  colnames(table_loop_r) <- c("N","T","Mean β1=1","SD β1","Mean β2=3","SD β2","Mean μ=5","SD μ",
                                 "Mean γ=2","SD γ","Mean =4","SD δ")
  table_loop_r_list[[i]] <- table_loop_r
}

