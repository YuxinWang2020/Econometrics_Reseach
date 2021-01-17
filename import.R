# arrayToList <- function(X){
#   X_list <- list()
#   T_ <- dim(X)[1]
#   N<-dim(X)[2]
#   k <- dim(X)[3]
#   for (i in 1:N){
#     X_list[[i]]<-X[,i,]
#   }
#   return(X_list)
# }

# write.csv(sim.data, file = "sim2.csv", row.names = FALSE)
beta <- read.csv("out/Archive(T100N100)/beta.txt", header = F, sep = "\t")
beta0 <- read.csv("out/Archive(T100N100)/beta0.txt", header = F, sep = "\t")
X_it1 <- read.csv("out/Archive(T100N100)/X_it1.txt", header = F, sep = "\t")   # nrow = T_
X_it2 <- read.csv("out/Archive(T100N100)/X_it2.txt", header = F, sep = "\t")
X_it3 <- read.csv("out/Archive(T100N100)/X_it3.txt", header = F, sep = "\t")
X_it4 <- read.csv("out/Archive(T100N100)/X_it4.txt", header = F, sep = "\t")
X_it5 <- read.csv("out/Archive(T100N100)/X_it5.txt", header = F, sep = "\t")
Y <- read.csv("out/Archive(T100N100)/Y.txt", header = F, sep = "\t")
# 
# X_array <- array(NA, dim = c(100,100,5))
# X_array[,,1] <- as.matrix(X_it1)
# X_array[,,2] <- as.matrix(X_it2)
# X_array[,,3] <- as.matrix(X_it3)
# X_array[,,4] <- as.matrix(X_it4)
# X_array[,,5] <- as.matrix(X_it5)
# Y_array <- as.matrix(Y)

T_ = 100
N = 100
p=5
DGP3 <- function(){
  # Set parameters
  T_ = 100
  N = 100
  p=5
  # Simulate data
  X_it_1 <- t(X_it1) # X_it_1  nrow=N, ncol=T_
  X_it_2 <- t(X_it2)
  X_it_3 <- t(X_it3)
  X_it_4 <- t(X_it4)
  X_it_5 <- t(X_it5)
  Y_it <- t(Y)
  # Save all results to data frame
  X_df <- data.frame(x_it_1 = as.vector(t(X_it_1)),
                     x_it_2 = as.vector(t(X_it_2)),
                     x_it_3 = as.vector(t(X_it_3)),
                     x_it_4 = as.vector(t(X_it_4)),
                     x_it_5 = as.vector(t(X_it_5)))
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   y_it = as.vector(t(Y_it)),
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
dgp3 <- DGP3()

least_squares <- function(X_list, Y_list, df, tolerance, beta_hat_0){
  # Initialize
  r <- 2
  K <- dim(X_list[[1]])[2]
  # formulate <- reformulate(response = c("y_it"), termlabels =
  #                            paste0("x_it_",c(1:K)) %>% paste(collapse = " + "))
  # beta_hat_0 <- plm(formulate, data=df, model="pooling")$coefficients %>% as.matrix()
  
  beta_hat_list <- list(beta_hat_0)
  e <- Inf
  while (e > tolerance) {
    F_hat <- calculate_F_hat(X_list, Y_list, beta_hat_0, r)
    Lambda_hat <- calculate_Lambda_hat(X_list, Y_list, beta_hat_0, F_hat, r)
    beta_hat <- calculate_beta_hat(X_list, Y_list, F_hat, Lambda_hat)
    beta_hat_list[[length(beta_hat_list)+1]] <- beta_hat
    e <- norm(beta_hat - beta_hat_0, type = "F")
    beta_hat_0 <- beta_hat
  }
  
  return(list(beta_hat=beta_hat, beta_hat_list=beta_hat_list, F_hat=F_hat, Lambda_hat=Lambda_hat))
}
ls <- least_squares(dgp3$X_list, dgp3$Y_list, dgp3$df, 0.0001, as.matrix(beta0))
beta
ls$beta_hat



