rm(list = ls())

if (!require("dplyr")) install.packages("dplyr")
if (!require("MASS")) install.packages("MASS")
if (!require("plm")) install.packages("plm")
if (!require("matlab")) install.packages("matlab") # for meshgrid()
# for plot
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("akima")) install.packages("akima")
if (!require("plotly")) install.packages("plotly")
# use JIT
library(compiler)
setCompilerOptions(optimize=3)
enableJIT(3)

# create dir
dir.create("../out", showWarnings = F)
dir.create("../out/figures", showWarnings = F)

set.seed(123)
source("DGPs.R")
source("Methods.R")
source("Statistics.R")

coefficients <- c("beta1", "beta2", "mu", "gamma", "delta")
select.coef <- 1 # choose 1~5 in coefficients

##### DGP1 #####

#### sim for different N ####
all_N <- c(10,20,30,40,50,100)
T_ <- 50
all_T <- rep(T_, length(all_N))
nsims <- 1000
beta_true <- c(1,3,5,2,4)
sim_figure_dgp1_1 <- sim_dgp1_fe(beta_true, all_N, all_T, nsims)
stat_figure_dgp1_1 <- statistics(sim_figure_dgp1_1$df_beta_hat_fe, NULL, beta_true, all_N, all_T, nsims)

### plot beta_hat for different N ###
select.col <- paste0("beta.",select.coef)
select.df_beta_hat <- select(sim_figure_dgp1_1$df_beta_hat_fe, c(1:3, select.col))
select.df_beta_hat$N <- as.factor(select.df_beta_hat$N)
plot.heigth <- max(beta_true[select.coef] - min(select.df_beta_hat[select.col]),
                   max(select.df_beta_hat[select.col]) - beta_true[select.coef])

## violin plot for beta_hat in range all_N
violin_plot <- ggplot(data = select.df_beta_hat, aes(x=N)) +
  geom_violin(aes(y = get(select.col), fill = N), color = "black", width=0.9, alpha=0.9) +
  geom_boxplot(aes(y = get(select.col)), alpha=0.6, fill=I("white"), width=0.06) +
  geom_hline(yintercept = beta_true[select.coef], color = I("black"), alpha=0.5) +
  ylim(beta_true[select.coef] - plot.heigth, beta_true[select.coef] + plot.heigth) +
  labs(title = paste0("fix T=",T_," (dgp1)"), y = paste(coefficients[select.coef], "hat")) +
  scale_fill_brewer(palette = "Blues") +
  guides(fill=FALSE, alpha=FALSE, color=FALSE, shape=FALSE) +
  theme_minimal()
print(violin_plot)
ggsave(paste0("../out/figures/dgp1_beta_hat_T",T_,"_violin.png"),width=5,height=5,units="in",dpi = 300)

#### sim for different N & T ####
all_N <- seq(10,100,10)
all_T <- seq(10,100,10)
grid <- meshgrid(all_N, all_T)
all_N <- as.vector(grid$x)
all_T <- as.vector(grid$y)

nsims <- 100
beta_true <- c(1,3,5,2,4)
sim_figure_dgp1_2 <- sim_dgp1_fe(beta_true, all_N, all_T, nsims)
stat_figure_dgp1_2 <- statistics(sim_figure_dgp1_2$df_beta_hat_fe, NULL, beta_true, all_N, all_T, nsims)

#### plot rmse for different N and T ####
select.col <- paste0("rmse.",select.coef)
# heat plot
im <- with(stat_figure_dgp1_2, interp(T_, N, get(select.col)))
with(im,image(x,y,z, xlab = "T_", ylab = "N"))
# point plot
point_plot_N <- ggplot(data = stat_figure_dgp1_2, aes(x=N)) +
  geom_jitter(aes(y = get(select.col), color = as.factor(N)), height=0, alpha=0.8) +
  geom_smooth(aes(y = get(select.col)), method = "loess", se=F, color=I("azure4"), formula = "y~x", size=0.5) +
  labs(title = paste0("fix T=",T_," (dgp1)"), x = "N", y = paste(coefficients[select.coef], "rmse")) +
  scale_color_brewer(palette = "Paired") +
  guides(fill=FALSE, alpha=FALSE, color=FALSE, shape=FALSE) +
  theme_minimal()
print(point_plot_N)
ggsave("../out/figures/dgp1_rmse_point.png",width=5,height=5,units="in",dpi = 300)



#### GDP2 ####

##### sim for different models & N for interactive-effect estimator and within estimator #####
all_N <- c(10,20,30,40,50,100)
T_ <- 100
all_T <- rep(T_, length(all_N))
nsims <- 1000
beta_true <- c(1,3,5,2,4)
tolerance <- 0.0001
r <- 2
models <- c("model1","model2","model3","model4")
sim_figure_dgp2_1 <- lapply(as.list(models),
                            function(model)
                              sim_dgp2_ls_fe(beta_true, tolerance, r, model, all_N, all_T, nsims, need.sde=F, need.fe=T))
names(sim_figure_dgp2_1) <- models
for(model in models){
  sim_data <- sim_figure_dgp2_1[[model]]
  ### plot beta_hat for different N ###
  select.col <- paste0("beta.",select.coef)
  select.df_beta_hat_ls <- select(sim_data$df_beta_hat_ls, c(1:3, select.col))
  select.df_beta_hat_fe <- select(sim_data$df_beta_hat_fe, c(1:3, select.col))
  select.df_beta_hat <- bind_rows(cbind(method="interactive-effect estimator", select.df_beta_hat_ls), cbind(method="within estimator", select.df_beta_hat_fe))
  select.df_beta_hat$N <- as.factor(select.df_beta_hat$N)
  select.df_beta_hat$method <- as.factor(select.df_beta_hat$method)
  
  # jitter point plot for beta_hat in range all_N
  point_plot <- ggplot(data = select.df_beta_hat, aes(x=N)) +
    geom_jitter(aes(y = get(select.col), color = method, shape = method), height=0) +
    geom_hline(yintercept = beta_true[select.coef], color = I("black"), alpha=0.5) +
    labs(title = paste0("fix T=",T_), x = "N", y = paste(coefficients[select.coef], "hat")) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal()
  print(point_plot)
  ggsave(paste0("../out/figures/dgp2_",model,"_beta_hat_T",T_,"_point.png"),width=7,height=5,units="in",dpi = 300)
  ## box plot for beta_hat in range all_N
  box_plot <- ggplot(data = select.df_beta_hat, aes(x=N)) +
    geom_boxplot(aes(y=get(select.col), fill = N), alpha=0.8, color=I("black"), width=0.6) +
    geom_hline(yintercept = beta_true[select.coef], color = I("black"), alpha=0.5) +
    labs(title = paste0("fix T=",T_), y = paste(coefficients[select.coef], "hat")) +
    scale_fill_brewer(palette = "Blues") +
    guides(fill=FALSE, alpha=FALSE, color=FALSE, shape=FALSE) +
    theme_minimal() +
    facet_wrap(~method)
  print(box_plot)
  ggsave(paste0("../out/figures/dgp2_",model,"_beta_hat_T",T_,"_box.png"),width=5,height=5,units="in",dpi = 300)
  ## violin plot for beta_hat in range all_N
  violin_plot <- ggplot(data = select.df_beta_hat, aes(x=N)) +
    geom_violin(aes(y = get(select.col), fill = N), color = "black", width=0.9, alpha=0.9) +
    geom_boxplot(aes(y = get(select.col)), alpha=0.6, fill=I("white"), width=0.06) +
    geom_hline(yintercept = beta_true[select.coef], color = I("black"), alpha=0.5) +
    labs(title = paste0("fix T=",T_), y = paste(coefficients[select.coef], "hat")) +
    scale_fill_brewer(palette = "Blues") +
    guides(fill=FALSE, alpha=FALSE, color=FALSE, shape=FALSE) +
    theme_minimal() +
    facet_wrap(~method)
  print(violin_plot)
  ggsave(paste0("../out/figures/dgp2_",model,"_beta_hat_T",T_,"_violin.png"),width=10,height=5,units="in",dpi = 300)
  
}



##### sim for different models & N & T for interactive-effect estimator and within estimator #####
all_N <- seq(10,100,10)
all_T <- seq(10,100,10)
grid <- meshgrid(all_N, all_T)
all_N_grid <- as.vector(grid$x)
all_T_grid <- as.vector(grid$y)
nsims <- 100
beta_true <- c(1,3,5,2,4)
tolerance <- 0.0001
r <- 2
models <- c("model1","model2","model3","model4")
sim_figure_dgp2_2 <- lapply(as.list(models),
                            function(model)
                              sim_dgp2_ls_fe(beta_true, tolerance, r, model, all_N_grid, all_T_grid, nsims, need.sde=F, need.fe=T))
stat_figure_dgp2_2_ls <- lapply(sim_figure_dgp2_2,
                                function(sim_data)
                                  statistics(sim_data$df_beta_hat_ls, sim_data$df_sde, beta_true, all_N_grid, all_T_grid, nsims))
stat_figure_dgp2_2_fe <- lapply(sim_figure_dgp2_2,
                                function(sim_data)
                                  statistics(sim_data$df_beta_hat_fe, NULL, beta_true, all_N_grid, all_T_grid, nsims))
names(stat_figure_dgp2_2_fe) <- names(stat_figure_dgp2_2_ls) <- names(sim_figure_dgp2_2) <- models

for(model in models){
  stat_ls <- stat_figure_dgp2_2_ls[[model]]
  stat_fe <- stat_figure_dgp2_2_fe[[model]]
  
  #### plot rmse for different N and T ####
  select.col <- paste0("rmse.",select.coef)
  
  # Heatmap ls (N, T, rmse)
  heatmap_ls <- ggplot(stat_ls, aes(x=N, y=T_, fill=get(select.col))) +
    geom_tile() +
    scale_fill_distiller(palette = "YlOrBr", direction = 1) +
    labs(title = paste0("interactive-effect estimator rmse"), x = "N", y = "T") +
    theme_minimal() +
    theme(legend.title = element_blank())
  print(heatmap_ls)
  ggsave(paste0("../out/figures/dgp2_",model,"_rmse_heatmap_ls.png"),width=5,height=5,units="in",dpi = 300)
  
  # Heatmap fe (N, T, rmse)
  heatmap_fe <- ggplot(stat_fe, aes(x=N, y=T_, fill=get(select.col))) +
    geom_tile() +
    scale_fill_distiller(palette = "YlOrBr", direction = 1) +
    labs(title = paste0("within estimator rmse"), x = "N", y = "T") +
    theme_minimal() +
    theme(legend.title = element_blank())
  print(heatmap_fe)
  ggsave(paste0("../out/figures/dgp2_",model,"_rmse_heatmap_fe.png"),width=5,height=5,units="in",dpi = 300)
  
  # point plot ls
  point_plot_N_ls <- ggplot(data = stat_ls, aes(x=N)) +
    geom_jitter(aes(y = get(select.col), color = as.factor(N)), height=0, alpha=0.8) +
    geom_smooth(aes(y = get(select.col)), method = "loess", se=F, color=I("azure4"), formula = "y~x", size=0.5) +
    labs(title = paste0("interactive-effect estimator"), x = "N", y = paste(coefficients[select.coef], "rmse")) +
    scale_color_brewer(palette = "Paired") +
    guides(fill=FALSE, alpha=FALSE, color=FALSE, shape=FALSE) +
    theme_minimal()
  print(point_plot_N_ls)
  ggsave(paste0("../out/figures/dgp2_",model,"_rmse_point_ls.png"),width=5,height=5,units="in",dpi = 300)
  
  # point plot fe
  point_plot_N_fe <- ggplot(data = stat_fe, aes(x=N)) +
    geom_jitter(aes(y = get(select.col), color = as.factor(N)), height=0, alpha=0.8) +
    geom_smooth(aes(y = get(select.col)), method = "loess", se=F, color=I("azure4"), formula = "y~x", size=0.5) +
    labs(title = paste0("within estimator"), x = "N", y = paste(coefficients[select.coef], "rmse")) +
    scale_color_brewer(palette = "Paired") +
    guides(fill=FALSE, alpha=FALSE, color=FALSE, shape=FALSE) +
    theme_minimal()
  print(point_plot_N_fe)
  ggsave(paste0("../out/figures/dgp2_",model,"_rmse_point_fe.png"),width=5,height=5,units="in",dpi = 300)
}

##### sim for different r for ls #####
r_N <- c(50)
r_T <- c(50)
nsims <- 1000
beta_true <- c(1,3,5,2,4)
tolerance <- 0.0001
model <- "model4"
rs <- c(1:10)
sim_figure_dgp2_3 <- as.list(rep(NA, length(rs)))
sim_figure_dgp2_3 <- lapply(as.list(rs),
                            function(r)
                              sim_dgp2_ls_fe(beta_true, tolerance, r, model, r_N, r_T, nsims, need.sde=F, need.fe=F))
stat_figure_dgp2_3_ls <- lapply(sim_figure_dgp2_3,
                                function(sim_data)
                                  statistics(sim_data$df_beta_hat_ls, sim_data$df_sde, beta_true, r_N, r_T, nsims))
names(stat_figure_dgp2_3_ls) <- names(sim_figure_dgp2_3) <- paste0("r",rs)

stat_figure_loop_r <- Reduce(function(df, r) bind_rows(df,  cbind(r=r, stat_figure_dgp2_3_ls[[paste0("r",r)]])), rs, init=data.frame(matrix(NA, nrow = 0, ncol = 0)))
beta_hat_figure_loop_r <- Reduce(function(df, r) bind_rows(df,  cbind(r=r, sim_figure_dgp2_3[[paste0("r",r)]]$df_beta_hat_ls)), rs, init=data.frame(matrix(NA, nrow = 0, ncol = 0)))


# point plot rmse
select.col <- paste0("rmse.",select.coef)
point_plot_r <- ggplot(data = filter(stat_figure_loop_r,r!=1), aes(x=r)) +
  geom_jitter(aes(y = get(select.col), color = as.factor(r)),height=0,width=0, alpha=0.8) +
  geom_smooth(aes(y = get(select.col)), method = "loess", se=F, color=I("azure4"), formula = "y~x", size=0.5) +
  labs(title = paste0("fix N=",r_N," T=",r_T," (dgp2 ",model,")"), x = "r", y = paste(coefficients[select.coef], "rmse")) +
  scale_color_brewer(palette = "Paired") +
  guides(fill=FALSE, alpha=FALSE, color=FALSE, shape=FALSE) +
  theme_minimal()
print(point_plot_r)
ggsave(paste0("../out/figures/dgp2_r_rmse_point.png"),width=5,height=5,units="in",dpi = 300)
# box plot for beta_hat
select.col <- paste0("beta.",select.coef)
box_plot <- ggplot(data = beta_hat_figure_loop_r, aes(x=as.factor(r))) +
  geom_boxplot(aes(y=get(select.col), fill = as.factor(r)), alpha=0.8, color=I("black"), width=0.6) +
  geom_hline(yintercept = beta_true[select.coef], color = I("black"), alpha=0.5) +
  labs(title = paste0("fix N=",r_N," T=",r_T," (dgp2 ",model,")"), y = paste(coefficients[select.coef], " hat"), x="r") +
  scale_fill_brewer(palette = "Spectral") +
  guides(fill=FALSE, alpha=FALSE, color=FALSE, shape=FALSE) +
  theme_minimal()
print(box_plot)
ggsave(paste0("../out/figures/dgp2_r_beta_box.png"),width=5,height=5,units="in",dpi = 300)


save.image(file = "../out/figures/figures.RData")
