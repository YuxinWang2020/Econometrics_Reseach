#######################
#####     DGP     #####
#######################


# Y_i <- alpha_i + X_i %*% beta_true + u_i
DGP1 <- function(T_, N, beta_true){
  # Set parameters
  p <- length(beta_true)
  range_x <- c(0,1)
  alpha_i <- 1
  mu_u_i <- 0
  sigma_u_i <- 1
  # Simulate data 
  X <- array(runif(p*T_*N,range_x[1],range_x[2]), dim = c(T_,N,p)) # array of c(T,N,p)
  U <- array(rnorm(N*T_, mean=mu_u_i, sd=sigma_u_i), dim = c(T_,N)) # matrix of c(T,N)
  Y <- Reduce(function(s,i) s + beta_true[i] * X[,,i], 1:length(beta_true), init=0) + U + alpha_i # matrix of c(T,N)
  X_list <- lapply(seq(dim(X)[2]), function(x) X[ , x, ])
  Y_list <- lapply(seq(dim(Y)[2]), function(x) Y[ , x])
  # Save all results to data frame
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   y_it = as.vector(Y),
                   x_it = matrix(X, nrow=T_*N, ncol=p))
  colnames(df)[4:(3+p)] <- paste0("x_it_",c(1:(p)))
  
  return(list(X_list = X_list, Y_list = Y_list, df=df))
}

# model5: y_it = mu + beta1*x_it_1 + beta2*x_it_2 + x_i*gamma + w_t*delta + Lambda_i*Factor_t + eps
# model4: y_it = mu + beta1*x_it_1 + beta2*x_it_2 + Lambda_i*Factor_t + eps
# model3: y_it = mu + beta1*x_it_1 + beta2*x_it_2 + alpha_i + z_t + eps
# model2: y_it = mu + beta1*x_it_1 + beta2*x_it_2 + alpha_i + eps
# model1: y_it = beta1*x_it_1 + beta2*x_it_2 + alpha_i + eps

# beta_true = c(beta1, beta2, mu, gamma, delta)
DGP2 <- function(T_, N, beta_true, model){
  p <- ifelse(model == "model5", 5, ifelse(model == "model1", 2, 3))
  # Set parameters
  if(model %in% c("model5", "model4", "model3", "model2")){
    mu <- beta_true[3]
  } else{
    mu <- 0
  }
  if(model == "model5"){
    gamma <- beta_true[4]
    delta <- beta_true[5]
  } else {
    gamma <- delta <- 0
  }
  iota <- rep(1, 2)
  mu1 <- mu2 <- c1 <- c2 <- 1
  # Generate variables
  if(model %in% c("model5", "model4")){
    Factor <- matrix(rnorm(n = 2*T_, mean = 0, sd = 1), nrow=2, ncol=T_)
    Lambda <- matrix(rnorm(n = 2*N, mean = 0, sd = 1), nrow=2, ncol=N)
  } else if(model == "model3"){
    Factor <- rbind(matrix(rnorm(n = T_, mean = 0, sd = 1), nrow=1, ncol=T_), 1)
    Lambda <- rbind(1, matrix(rnorm(n = N, mean = 0, sd = 1), nrow=1, ncol=N))
  } else{
    Factor <- rbind(0, matrix(1, nrow=1, ncol=T_))
    Lambda <- rbind(1, matrix(rnorm(n = N, mean = 0, sd = 1), nrow=1, ncol=N))
  }
  Eta_1 <- matrix(rnorm(n = T_*N, mean = 0, sd = 1), nrow=N, ncol=T_)
  Eta_2 <- matrix(rnorm(n = T_*N, mean = 0, sd = 1), nrow=N, ncol=T_)
  Eps <- matrix(rnorm(n = T_*N, mean = 0, sd = 2), nrow=N, ncol=T_)
  e <- rnorm(n = N, mean = 0, sd = 1)
  eta <- rnorm(n = T_, mean = 0, sd = 1)
  # Calculate intermediate variables
  iota_Lambda <- crossprod(Lambda, iota)
  iota_Factor <- crossprod(iota, Factor)
  Lambda_Factor <- crossprod(Lambda, Factor)
  x <- iota_Lambda + e
  w <- iota_Factor + eta
  # Simulate data
  X_1 <- mu1 + c1 * Lambda_Factor + iota_Lambda %*% rep(1, T_) + rep(1, N) %*% iota_Factor + Eta_1
  X_2 <- mu2 + c2 * Lambda_Factor + iota_Lambda %*% rep(1, T_) + rep(1, N) %*% iota_Factor + Eta_2
  X_3 <- matrix(1, nrow=N, ncol=T_)
  X_4 <- x %*% rep(1, T_)
  X_5 <- rep(1, N) %*% w
  Y <- beta_true[1]*X_1 + beta_true[2]*X_2 + mu*X_3 + gamma*X_4 + delta*X_5 + Lambda_Factor + Eps
  # Save all results to data frame
  X_df <- data.frame(x_it_1 = as.vector(t(X_1)),
                     x_it_2 = as.vector(t(X_2)),
                     x_it_3 = as.vector(t(X_3)),
                     x_it_4 = as.vector(t(X_4)),
                     x_it_5 = as.vector(t(X_5)))
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   y_it = as.vector(t(Y)),
                   X_df[,1:p])
  # Save results to lists
  X_list <- list()
  Y_list <- list()
  for(i in 1:N){
    X_list[[i]] <- as.matrix(df[((i-1)*T_+1):(i*T_), 4:(3+p)])
    Y_list[[i]] <- as.matrix(df$y_it[((i-1)*T_+1):(i*T_)])
  }
  return(list(df=df, X_list=X_list, Y_list=Y_list))
}
