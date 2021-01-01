rm(list = ls())
if (!require("MASS")) install.packages("MASS")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("akima")) install.packages("akima")
if (!require("plotly")) install.packages("plotly")
if (!require("reshape2")) install.packages("reshape2")
# if (!require("cowplot")) install.packages("cowplot")
if (!require("plm")) install.packages("plm")
set.seed(100)


##### DGP #####
DGP <- function(K, T_, N, beta_i){
  if(length(beta_i) != K){
    warning("length of beta_i != K")
    return(c())
  }
  range_x = c(0,1)
  alpha_i = 1
  mu_u_i = rep(0, T_)
  sigma_u_i = diag(rep(1, T_))
  
  X <- lapply(as.list(1:N), function(i) array(runif(K*T_,range_x[1],range_x[2]), dim = c(T_, K)))
  U <- lapply(as.list(1:N), function(i) mvrnorm(n = 1, mu = mu_u_i, Sigma = sigma_u_i))
  Y <- mapply(function(X_i, u_i) alpha_i + X_i %*% beta_i + u_i, X, U, SIMPLIFY=F)
  # X <- list()
  # Y <- list()
  # for(i in 1:N){
  #     X_i <- array(runif(K*T_,range_x[1],range_x[2]), dim = c(T_, K))
  #     u_i = mvrnorm(n = 1, mu = mu_u_i, Sigma = sigma_u_i)
  #     Y_i = alpha_i + X_i %*% beta_i + u_i
  #     
  #     X[[i]] = X_i
  #     Y[[i]] = Y_i
  # }
  
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   y_it = NA) %>% cbind(matrix(NA,nrow=N*T_,ncol=K))
  colnames(df)[4:(3+K)] <- paste("x",c(1:K), sep = "")

  for(i in 1:N){
    for(t in 1:T_){
      df[(i-1)*T_ + t, 4:(3+K)] <- X[[i]][t,]
      df$y_it[(i-1)*T_ + t] <- Y[[i]][t]
    }
  }
  return(list(X = X, Y = Y, df=df))
}

DGP2 <- function(T_, N, beta=c(1,2)){
  if(T_ < 2 || N < 2)
    stop("T or N less than 2")
  
  alpha_i <- runif(N, min = 0, max = 1)
  eps_t <- rnorm(T_, mean = 0, sd = 1)
  z0 <- 0
  z_t <- c(0.5*z0 + eps_t[1])
  for(i in 2:T_)
    z_t[i] <- 0.5*z_t[i-1] + eps_t[i]
  v_i <- rnorm(N, mean = 0, sd = 2)
  e_t <- rnorm(T_ + 1, mean = 0, sd = 1)
  mu <- 5
  w_t <- mu + e_t[2:(T_+1)] + e_t[1:T_] / 3
  x0 = 1
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   x0 = x0,
                   x_it = NA,
                   u_it = NA,
                   y_it = NA)
  for(i in 1:N){
    for(t in 1:T_){
      df$x_it[(i-1)*T_ + t] <- sqrt(alpha_i[i]) * z_t[t]
      df$u_it[(i-1)*T_ + t] <- v_i[i] + w_t[t]
      df$y_it[(i-1)*T_ + t] <- alpha_i[i] + beta[1] * x0 +
        beta[2] * df$x_it[(i-1)*T_ + t] + df$u_it[(i-1)*T_ + t]
    }
  }
  return(df)
}

OLS_FE <- function(X, Y){
  if(length(dim(X[[1]])) != 2){
    warning("X_i is not 2d array")
    return(c())
  }
  
  N = length(X)
  K = dim(X[[1]])[2]
  T_ = dim(X[[1]])[1]
  A = array(0, dim = c(K,K))
  B = rep(0, K)
  for(i in 1:N){
    x_i_mean = colMeans(X[[i]])
    y_i_mean = mean(Y[[i]])
    for(t in 1:T_){
      A = A + (X[[i]][t,] - x_i_mean) %*% t(X[[i]][t,] - x_i_mean)
      B = B + (X[[i]][t,] - x_i_mean) * (Y[[i]][t] - y_i_mean)
    }
  }
  beta_hat_fe = solve(A) %*% B
  return(list(beta_hat = beta_hat_fe))
}

OLS_FE2 <- function(df){
  data <- pdata.frame(df,index=c("i","t"))
  result <- plm(y_it ~ x0 + x_it, data=data,
                effect = "twoways",model="within")
  return(result$coefficients)
}

##### test #####
df <- DGP2(5, 7)
OLS_FE2(df)

dgp = DGP(K = 3,
          T_ = 3,
          N = 400,
          beta_i = c(1,2,3))
OLS_FE(dgp$X, dgp$Y)
plm(y_it ~ x1+x2+x3, data=dgp$df,
    effect = "twoways",model="within")$coefficients


##### simulation #####

sim <- function(K, beta_true, all_n, all_T, nsims){
  beta_hat_1 = array(NA, dim = c(length(all_n), length(all_T), nsims))
  summary_beta_hat_1 = data.frame(T_ = rep(all_T, times = length(all_n)),
                               N = rep(all_n, each = length(all_T)),
                               mean = NA,
                               var = NA,
                               bias = NA)
  for(i in 1:length(all_n)){
    N = all_n[i]
    for(j in 1:length(all_T)){
      T_ = all_T[j]
      for(h in 1:nsims){
        dgp = DGP(K = K,
                  T_ =T_,
                  N = N,
                  beta_i = beta_true)
        result = OLS_FE(dgp$X, dgp$Y)
        beta_hat = result$beta_hat
        beta_hat_1[i, j, h] = beta_hat[1]
      }
      summary_beta_hat_1$mean[(i-1)*length(all_T)+j] = mean(beta_hat_1[i, j,])
      summary_beta_hat_1$var[(i-1)*length(all_T)+j] = var(beta_hat_1[i, j,])
      summary_beta_hat_1$bias[(i-1)*length(all_T)+j] =
        abs(summary_beta_hat_1$mean[(i-1)*length(all_T)+j] - beta_true[1]) / beta_true[1]
    }
  }
  return(list(beta_hat_1=beta_hat_1, summary_beta_hat_1=summary_beta_hat_1))
}

##### sim1 #####
all_n = c(1000, 2000, 3000)
nsims=100
beta_true = c(1:3)
sim1 <- sim(K = 3, beta_true = beta_true, all_n = all_n, all_T = c(3), nsims = nsims)
beta_hat_1 = sim1$beta_hat_1
beta_true_1 = beta_true[1]
beta_hat_1_df = data.frame(t(beta_hat_1[,1,]))
colnames(beta_hat_1_df) <- all_n
beta_hat_1_df <- melt(beta_hat_1_df, measure.vars=colnames(beta_hat_1_df), variable.name = "N")
# generate point and box plot for every N
plot.heigth <- max(beta_true_1 - min(beta_hat_1_df$value), max(beta_hat_1_df$value) - beta_true_1)
point_plot <- ggplot(data = beta_hat_1_df) +
  geom_point(aes(y = value, x = c(1:(nsims*length(all_n))), color = N, shape = N)) +
  geom_hline(yintercept = beta_true_1, color = I("black")) +
  ylim(beta_true_1 - plot.heigth, beta_true_1 + plot.heigth) +
  labs(title = "point plot for every N", x = "iteration", y = "beta_hat_1") +
  theme(axis.text.x = element_blank())
print(point_plot)

box_plot <- ggplot(data = beta_hat_1_df) +
  geom_boxplot(aes(y = value, x = N, color = N)) +
  ylim(min(beta_hat_1_df$value), max(beta_hat_1_df$value)) +
  labs(title = "box plot for every N", y = "beta_hat_1")
print(box_plot)

##### sim2 #####
sim2 <- sim(K = 3, beta_true = c(1:3), all_n = seq(1000, 10000, 1000), all_T = seq(3,30,3), nsims = 1)
beta_hat_1 = sim2$beta_hat_1
beta_true_1 = 1
summary_beta_hat_1 = sim2$summary_beta_hat_1

im <- with(summary_beta_hat_1, interp(T_,N,bias))
with(im,image(x,y,z, xlab = "T_", ylab = "N"))

plot_ly(x=im$x, y=im$y, z=im$z, type="surface") %>% add_surface() %>% layout(scene = list(
    xaxis = list(title = 'T'),
    yaxis = list(title = 'N'),
    zaxis = list(title = 'bias')))
