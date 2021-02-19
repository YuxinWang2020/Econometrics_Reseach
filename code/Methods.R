#####################
####   Method    ####
#####################

#####  Simple Fixed Effects Model   #####

# For model 1, we write a function to calculate the within estimator.
OLS_FE <- function(X_list, Y_list){ 
  # Initialize
  N <- length(X_list)
  T_ <- dim(X_list[[1]])[1]
  p <- dim(X_list[[1]])[2]
  A <- array(0, dim = c(p,p))
  B <- rep(0, p)
  # Loop over N & T_
  for(i in 1:N){
    X_i <- X_list[[i]]
    Y_i <- Y_list[[i]]
    x_i_mean <- colMeans(X_i)
    y_i_mean <- mean(Y_i)
    
    A <- A + tcrossprod(t(X_i)-x_i_mean)
    B <- B + (t(X_i)-x_i_mean) %*% (Y_i - y_i_mean)
  }
  beta_hat_fe <- solve(A) %*% B
  return(list(beta_hat = beta_hat_fe))
}

# Alternatively, we can use plm package for estimation.
OLS_FE2 <- function(df){
  p <- ncol(df) - 3
  data <- pdata.frame(df,index=c("i","t"))
  formulate <- paste0("y_it ~ ", paste(paste0("x_it_",c(1:p)), collapse = " + "))
  result <- plm(formulate, data=data,
                effect = "individual",model="within")
  return(list(beta_hat = as.matrix(result$coefficients)))
}

# Now we use plm package to calculate model 2-4.
OLS_FE3 <- function(df){
  p <- ncol(df) - 3
  formulate <- paste0("y_it ~ ", paste(paste0("x_it_",c(1:p)), collapse = " + "))
  result <- plm(formulate, data=df,index=c("i","t"),
                effect = "twoways",model="within")
  return(list(beta_hat = as.matrix(result$coefficients)))
}



#####  Interacctive Fixed Effect Methods  #####

#Step 1:define funtion to calculate F_hat, dim of F_hat is (T_, r).
calculate_F_hat <- function(X_list, Y_list, beta_hat, r){
  N <- length(X_list)
  T_ <- dim(X_list[[1]])[1]
  p <- dim(X_list[[1]])[2]
  
  WWT <- matrix(0, nrow=T_, ncol=T_)
  for(i in 1:N){
    X_i <- X_list[[i]]
    Y_i <- Y_list[[i]]
    W_i <- (Y_i - X_i %*% beta_hat)
    WWT <- WWT + W_i %*% t(W_i)
  }
  eig <- eigen(WWT)
  F_hat <- as.matrix(sqrt(T_)*eig$vectors[,order(eig$values, decreasing = T)[1:r]])
  return(F_hat)
}

#Step 2:define funtion to calculate Lambda_hat, dim of Lambda_hat is (r, N).
calculate_Lambda_hat <- function(X_list, Y_list, beta_hat, F_hat, r){
  N <- length(X_list)
  T_ <- dim(X_list[[1]])[1]
  p <- dim(X_list[[1]])[2]
  
  Lambda_hat <- matrix(NA, nrow = N, ncol = r)
  for(i in 1:N){
    X_i <- X_list[[i]]
    Y_i <- Y_list[[i]]
    Lambda_hat[i,] <- t(F_hat) %*% (Y_i - X_i %*% beta_hat) / T_
  }
  return(Lambda_hat)
}

#Step 3:define funtion to calculate Beta_hat.
calculate_beta_hat <- function(X_list, Y_list, F_, Lambda){
  N <- length(X_list)
  T_ <- dim(X_list[[1]])[1]
  p <- dim(X_list[[1]])[2]
  
  A <- matrix(0, nrow=p, ncol=p)
  B <- matrix(0, nrow=p, ncol=1)
  for(i in 1:N){
    X_i <- X_list[[i]]
    Y_i <- Y_list[[i]]
    A <- A + t(X_i) %*% X_i
    B <- B + t(X_i) %*% (Y_i - F_ %*% Lambda[i,])
  }
  beta_hat <- solve(A) %*% B
  return(beta_hat)
}

#Step 4:calculate Beta_hat by iterations.
least_squares <- function(X_list, Y_list, df, tolerance, r, model){
  # Initialize
  p <- dim(X_list[[1]])[2]
  formulate <- paste0("y_it ~ ", paste0("x_it_",c(1:p)) %>% paste(collapse = " + "), ifelse(p<=2, " + 0", ""))
  if(model == "model2"){
    beta_hat_0 <- plm(formulate, data=df,effect = "twoways", model="within")$coefficients %>% as.matrix()
  } else{
    beta_hat_0 <- plm(formulate, data=df, model="pooling")$coefficients %>% as.matrix()
  }
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
