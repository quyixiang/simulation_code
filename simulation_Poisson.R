library(MASS)
#library(matlib)

simulate_poission <- function(p, n, m, seed) {

  # generate X[1,m+n], Y[1,n]
  X_n <- matrix(ncol = p, nrow = 0)
  set.seed(seed)
  X_add <- matrix(rpois(m*p, 10), ncol = p)
  Y_n <- c()
  set.seed(seed)
  for (i in c(1:(n))) {
    X_tmp <- matrix(rpois(p, 10), nrow=1)

      X_n <- rbind(X_n, X_tmp)
      Y_tmp <- rpois(1, 10*X_tmp[1])
      Y_n <- c(Y_n, Y_tmp)
  }
  
  # calculate beta_hat
  X_vec <- cbind(rep(1, nrow(X_n)), X_n)
  Y <- matrix(Y_n, nrow = length(Y_n))
  beta_hat <- solve(t(X_vec) %*% X_vec) %*% t(X_vec) %*% Y
  beta_2_hat <- beta_hat[c(2:nrow(beta_hat)), 1]
  
  # calculate Y_bar and X_bar
  Y_bar <- mean(Y_n)
  X_bar_n <- apply(X_n, 2, mean)
  
  X_all <- rbind(X_n, X_add)
  X_bar_all <- apply(X_all, 2, mean)
  
  # calculate theta_hat_LS
  theta_hat_LS <- Y_bar - sum(beta_2_hat * (X_bar_n-10))
  theta_hat_SSLS <- Y_bar - sum(beta_2_hat * (X_bar_n - X_bar_all))
  
  # calculate MSE and sigma_2_Y
  MSE <- sum((Y - X_vec %*% beta_hat)^2) / (n - p - 1)
  sigma_2_Y <- sum((Y - Y_bar)^2) / (n - 1)


  returnlist <- list("theta_hat_LS" = theta_hat_LS, "theta_hat_SSLS" = theta_hat_SSLS, "Y_bar" = Y_bar, "theta" = 100, "MSE" = MSE, "sigma_2_Y" = sigma_2_Y)
  return(returnlist)
}


results <- data.frame(matrix(nrow = 0, ncol = 9))
colnames(results) <- c("p", "m", "n", "l2_Y_bar", "l2_theta_hat_SSLS", "l2_theta_hat_LS", "LCI_Y_bar", "LCI_theta_hat_SSLS", "LCI_theta_hat_LS")

pnm_list <- list(c(1, 100, 100), c(1, 100, 1000), c(10, 100, 100), c(10, 100, 1000), c(50, 100, 100), c(50, 100, 1000), c(10, 500, 100), c(10, 500, 1000), c(50, 500, 100), c(50, 500, 1000), c(200, 500, 100), c(200, 500, 1000))

for (pnm in pnm_list) {
  p <- pnm[1]
  n <- pnm[2]
  m <- pnm[3]
  print(paste0("p:", p, ",n:", n, ",m:", m))
  l2_Y_bar <- c()
  l2_theta_hat_SSLS <- c()
  l2_theta_hat_LS <- c()
  LCI_Y_bar <- c()
  LCI_theta_hat_SSLS <- c()
  LCI_theta_hat_LS <- c()
  for (i in c(1:500)) {
    if (i %% 100 == 0) print(i)
    
    seed <- i^2
    tmp_simulation <- simulate_poission(p, n, m, seed)
    tmp_l2_Y_bar <- (tmp_simulation$Y_bar - tmp_simulation$theta)^2
    tmp_l2_theta_hat_SSLS <- (tmp_simulation$theta_hat_SSLS - tmp_simulation$theta)^2
    tmp_l2_theta_hat_LS <- (tmp_simulation$theta_hat_LS - tmp_simulation$theta)^2
    l2_Y_bar <- c(l2_Y_bar, tmp_l2_Y_bar)
    l2_theta_hat_SSLS <- c(l2_theta_hat_SSLS, tmp_l2_theta_hat_SSLS)
    l2_theta_hat_LS <- c(l2_theta_hat_LS, tmp_l2_theta_hat_LS)

    Z <- qnorm(0.975)
    tmp_LCI_Y_bar <- 2 * Z * sqrt(tmp_simulation$sigma_2_Y / n)
    tmp_LCI_theta_hat_SSLS <- 2 * Z * sqrt((m * tmp_simulation$MSE + n * tmp_simulation$sigma_2_Y) / (n * (m + n)))
    tmp_LCI_theta_hat_LS <- 2 * Z * sqrt(tmp_simulation$MSE / n)
    LCI_Y_bar <- c(LCI_Y_bar, tmp_LCI_Y_bar)
    LCI_theta_hat_SSLS <- c(LCI_theta_hat_SSLS, tmp_LCI_theta_hat_SSLS)
    LCI_theta_hat_LS <- c(LCI_theta_hat_LS, tmp_LCI_theta_hat_LS)
  }
  results[nrow(results) + 1, ] <- c(p, m, n, mean(l2_Y_bar), mean(l2_theta_hat_SSLS), mean(l2_theta_hat_LS), mean(LCI_Y_bar), mean(LCI_theta_hat_SSLS), mean(LCI_theta_hat_LS))
}

write.csv(results, "results_poisson.csv")
