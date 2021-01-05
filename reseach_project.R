rm(list = ls())

if (!require("MASS")) install.packages("MASS")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("akima")) install.packages("akima")
if (!require("plotly")) install.packages("plotly")
if (!require("reshape2")) install.packages("reshape2")
if (!require("plm")) install.packages("plm")

set.seed(100)

#######################
#####     DGP     #####
#######################

DGP1 <- function(K, T_, N, beta){
  if(length(beta) != K){
    warning("length of beta != K")
    return(c())
  }
  range_x = c(0,1)
  alpha_i = 1
  mu_u_i = rep(0, T_)
  sigma_u_i = diag(rep(1, T_))
  
  X <- lapply(as.list(1:N), function(i) array(runif(K*T_,range_x[1],range_x[2]), dim = c(T_, K)))
  U <- lapply(as.list(1:N), function(i) mvrnorm(n = 1, mu = mu_u_i, Sigma = sigma_u_i))
  Y <- mapply(function(X_i, u_i) alpha_i + X_i %*% beta + u_i, X, U, SIMPLIFY=F)
  
  #     X_i <- array(runif(K*T_,range_x[1],range_x[2]), dim = c(T_, K))
  #     u_i = mvrnorm(n = 1, mu = mu_u_i, Sigma = sigma_u_i)
  #     Y_i = alpha_i + X_i %*% beta + u_i
  
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   y_it = NA) %>% cbind(matrix(NA,nrow=N*T_,ncol=K))
  colnames(df)[4:(3+K)] <- paste0("x_it_",c(1:K))
  
  for(i in 1:N){
    X_i <- X[[i]]
    u_i = U[[i]]
    Y_i = Y[[i]]
    for(t in 1:T_){
      df[(i-1)*T_ + t, 4:(3+K)] <- X_i[t,]
      df$y_it[(i-1)*T_ + t] <- Y_i[t]
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
                   y_it = NA,
                   x_it_1 = x0,
                   x_it_2 = NA)
  for(i in 1:N){
    for(t in 1:T_){
      x_it <- sqrt(alpha_i[i]) * z_t[t]
      df$x_it_2[(i-1)*T_ + t] <- x_it
      u_it <- v_i[i] + w_t[t]
      df$y_it[(i-1)*T_ + t] <- alpha_i[i] + beta[1] * x0 + beta[2] * x_it + u_it
    }
  }
  return(list(df=df))
}

DGP3 <- function(T_, N, beta=c(1,3)){
  r<-2
  K<-2
  mu <- 5
  gamma <- 2
  delta <- 4
  iota <- rep(1, r)
  mu1 <- mu2 <- c1 <- c2 <- 1
  
  F_t <- matrix(rnorm(n = r*T_, mean = 0, sd = 1), nrow=r, ncol=T_)
  Lambda_i <- matrix(rnorm(n = r*N, mean = 0, sd = 1), nrow=r, ncol=N)
  Eta_it_1 <- matrix(rnorm(n = T_*N, mean = 0, sd = 1), nrow=N, ncol=T_)
  Eta_it_2 <- matrix(rnorm(n = T_*N, mean = 0, sd = 1), nrow=N, ncol=T_)
  Eps_it <- matrix(rnorm(n = T_*N, mean = 0, sd = 4), nrow=N, ncol=T_)
  e_i <- rnorm(n = N, mean = 0, sd = 1)
  eta_t <- rnorm(n = T_, mean = 0, sd = 1)
  
  x_i <- t(iota) %*% Lambda_i + e_i
  w_t <- t(iota) %*% F_t + eta_t
  
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   y_it = NA, x_it_1 = NA, x_it_2 = NA)
  X_list <- list()
  Y_list <- list()
  for(i in 1:N){
    for(t in 1:T_){
      X_it_1 <- mu1 + c1 * t(Lambda_i[,i]) %*% F_t[,t] + 
        t(iota) %*% Lambda_i[,i] + t(iota) %*% F_t[,t] + Eta_it_1[i,t]
      X_it_2 <- mu2 + c2 * t(Lambda_i[,i]) %*% F_t[,t] + 
        t(iota) %*% Lambda_i[,i] + t(iota) %*% F_t[,t] + Eta_it_2[i,t]
      Y_it <- X_it_1*beta[1] + X_it_2*beta[2] + mu + 
        x_i[i]*gamma + w_t[t]*delta + t(Lambda_i[,i])%*%F_t[,t] + Eps_it[i,t]
      
      df$x_it_1[(i-1)*T_ + t] <- X_it_1
      df$x_it_2[(i-1)*T_ + t] <- X_it_2
      df$y_it[(i-1)*T_ + t] <- Y_it
    }
    X_list[[i]] <- as.matrix(df[((i-1)*T_+1):(i*T_), 4:(3+K)])
    Y_list[[i]] <- as.matrix(df$y_it[((i-1)*T_+1):(i*T_)])
  }
  
  return(list(df=df, F_=F_t, Lambda=Lambda_i,
              Eta_it_1=Eta_it_1, Eta_it_2=Eta_it_2,
              Eps_it=Eps_it, e_i=e_i, eta_t=eta_t,
              X_list=X_list, Y_list=Y_list))
}

#####################
####   Method    ####
#####################

OLS_FE <- function(df, K){
  N = length(unique(df$i))
  T_ = length(unique(df$t))
  A = array(0, dim = c(K,K))
  B = rep(0, K)
  for(i in 1:N){
    X_i <- as.matrix(df[((i-1)*T_+1):(i*T_), 4:(3+K)])
    Y_i <- as.matrix(df$y_it[((i-1)*T_+1):(i*T_)])
    x_i_mean = colMeans(X_i)
    y_i_mean = mean(Y_i)
    for(t in 1:T_){
      A = A + (X_i[t,] - x_i_mean) %*% t(X_i[t,] - x_i_mean)
      B = B + (X_i[t,] - x_i_mean) * (Y_i[t] - y_i_mean)
    }
  }
  beta_hat_fe = solve(A) %*% B
  return(list(beta_hat = beta_hat_fe))
}

OLS_FE2 <- function(df, K){
  data <- pdata.frame(df,index=c("i","t"))
  formulate <- reformulate(response = c("y_it"), termlabels =
                             paste0("x_it_",c(1:K)) %>%
                             paste0("I(", ., ")") %>%
                             paste(collapse = " + "))
  result <- plm(formulate, data=data,
                effect = "individual",model="within")
  return(result$coefficients)
}

least_squares <- function(df, K=2, r=2, tolerance=0.0001){
  # Initialize
  N = length(unique(df$i))
  T_ = length(unique(df$t))
  # step1: define F0
  F_0 <- diag(sqrt(T_), nrow = T_, ncol = r)
  # step2: caculate Beta_hat(F)
  caculate_beta_hat_F <- function(F_){
    M_F <- diag(1, nrow = T_) - (F_ %*% t(F_)) / T_
    A <- matrix(0, nrow=K, ncol=K)
    B <- matrix(0, nrow=K, ncol=1)
    for(i in 1:N){
      X_i <- as.matrix(df[((i-1)*T_+1):(i*T_), 4:(3+K)])
      Y_i <- as.matrix(df$y_it[((i-1)*T_+1):(i*T_)])
      A <- A + t(X_i) %*% M_F %*% X_i
      B <- B + t(X_i) %*% M_F %*% Y_i
    }
    beta_hat_F <- solve(A) %*% B
    return(beta_hat_F)
  }
  # step3: caculate F_hat
  caculate_F_hat <- function(beta_hat_F){
    WWT <- matrix(0, nrow=T_, ncol=T_)
    for(i in 1:N){
      X_i <- as.matrix(df[((i-1)*T_+1):(i*T_), 4:(3+K)])
      Y_i <- as.matrix(df$y_it[((i-1)*T_+1):(i*T_)])
      W_i <- (Y_i - X_i %*% beta_hat_F)
      WWT <- WWT + W_i %*% t(W_i)
    }
    eig <- eigen(WWT)
    F_hat <- sqrt(T_)*eig$vectors[,order(eig$values, decreasing = T)[1:r]]
    return(F_hat)
  }
  # step4: caculate Beta_hat(F_hat) by iterations
  e <- Inf
  beta_hat_F0 <- caculate_beta_hat_F(F_0)
  beta_hat_list <- list(beta_hat_F0)
  while (e > tolerance) {
    F_hat <- caculate_F_hat(beta_hat_F0)
    beta_hat_F_hat <- caculate_beta_hat_F(F_hat)
    
    beta_hat_list[[length(beta_hat_list)+1]] <- beta_hat_F_hat
    e <- norm(beta_hat_F_hat - beta_hat_F0, type = "F")
    beta_hat_F0 <- beta_hat_F_hat
  }
  return(list(beta_hat=beta_hat_F_hat, beta_hat_list=beta_hat_list))
}


##### test #####
dgp2 <- DGP2(10, 100)
OLS_FE2(dgp2$df, K=2)

dgp1 = DGP1(K = 3,
            T_ = 3,
            N = 400,  beta = c(1,2,3))
OLS_FE(dgp1$df, K=3)
OLS_FE2(dgp1$df, K=3)
T_ <- 10
N <- 100
dgp3 <- DGP3(T_,N,beta=c(1,3))
df <-  dgp3$df
X_list <- dgp3$X_list
Y_list <- dgp3$Y_list

least_squares(df,r=2)$beta_hat
# plm(y_it ~ x_it_1+x_it_2, data=df,
#     effect = "twoways",model="within")$coefficients


###########################
####   Visualization   ####
###########################

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
        dgp = DGP1(K = K,
                  T_ =T_,
                  N = N,
                  beta = beta_true)
        result = OLS_FE(dgp$df,K)
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
