library(MASS)
library(randomForest)

#Bonferroni based method

p_X <- 20
sd_X <- 3
sd_Y <- 1
n_train <- 500

set.seed(12)
beta_1 <- runif(p_X)
beta_2 <- runif(p_X)
mu_X <- runif(p_X)

X_train <- mvrnorm(n_train, mu=mu_X, Sigma=sd_X*diag(p_X))
Y_train <- sapply(1:n_train, function(i) rnorm(1, mean=sum(beta_1*X_train[i,])+2*sin(sum(beta_2*X_train[i,])), sd=sd_Y*sqrt(sum(X_train[i,]^2))/p_X))
rf <- randomForest(X_train, Y_train)

n <- 1000
m <- 5
tn <- 500

alpha_set <- seq(0.05, 0.2, 0.025)
fd <- matrix(0, nrow=tn, ncol=length(alpha_set))
tp <- matrix(0, nrow=tn, ncol=length(alpha_set))
coverage <- matrix(0, nrow=tn, ncol=length(alpha_set))
width <- matrix(0, nrow=tn, ncol=length(alpha_set))

pb <- txtProgressBar(min=0, max=tn, style=3)

for(k in 1:tn) {
  X <- mvrnorm(n, mu=mu_X, Sigma=sd_X*diag(p_X))
  Y <- sapply(1:n, function(i) rnorm(1, mean=sum(beta_1*X[i,])+2*sin(sum(beta_2*X[i,])), sd=sd_Y*sqrt(sum(X[i,]^2))/p_X))
  muhat <- predict(rf, X)
  
  X_test <- mvrnorm(m, mu=mu_X, Sigma=sd_X*diag(p_X))
  Y_test <- sapply(1:m, function(i) rnorm(1, mean=sum(beta_1*X_test[i,])+2*sin(sum(beta_2*X_test[i,])), sd=sd_Y*sqrt(sum(X_test[i,]^2))/p_X))
  muhat_test <- predict(rf, X_test)
  
  n_1 <- floor(n/2)
  cal_abs_diff <- abs(muhat[1:n_1]-muhat[(n_1+1):(2*n_1)])
  cal_Y_null <- ((muhat[1:n_1]-muhat[(n_1+1):(2*n_1)])*(Y[1:n_1]-Y[(n_1+1):(2*n_1)]) <= 0)
  
  test_pairs <- t(combn(m, 2))
  test_abs_diff <- apply(test_pairs, 1, function(a) abs(muhat_test[a[1]]-muhat_test[a[2]]))
  test_null <- apply(test_pairs, 1, function(a) (muhat_test[a[1]]-muhat_test[a[2]])*(Y_test[a[1]]-Y_test[a[2]]) <= 0)
  test_non_null <- 1-test_null
  
  p_val <- sapply(test_abs_diff, function(v) (sum((v < cal_abs_diff)*cal_Y_null)+1)/(n_1+1))
  
  I_hat <- rank(muhat_test, ties.method="random")
  true_rank <- rank(Y_test[order(muhat_test)])
  
  for(j in seq_along(alpha_set)) {
    alpha <- alpha_set[j]
    rej_pairs <- p_val <= (2*alpha/(m*(m-1)))
    fd[k,j] <- sum(rej_pairs*test_null) > 0
    tp[k,j] <- sum(rej_pairs*test_non_null)/max(sum(test_non_null), 1)
    
    rank_interval <- cbind(rep(1, m), rep(m, m))
    rej_pairs_ind <- test_pairs[rej_pairs, , drop=FALSE]
    
    if(nrow(rej_pairs_ind) > 0) {
      rej_pairs_ind_ordered <- t(apply(rej_pairs_ind, 1, function(v) sort(I_hat[v])))
      for(t in 1:m) {
        lower_t <- sum(rej_pairs_ind_ordered[,2] == t)+1
        upper_t <- m-sum(rej_pairs_ind_ordered[,1] == t)
        rank_interval[t,] <- c(lower_t, upper_t)
      }
    }
    
    coverage[k,j] <- min((rank_interval[,1] <= true_rank)*(true_rank <= rank_interval[,2]))
    width[k,j] <- mean(rank_interval[,2]-rank_interval[,1]+1)
  }
  
  setTxtProgressBar(pb, k)
}
close(pb)

fwer <- colMeans(fd)
fwer_se <- apply(fd, 2, sd)/sqrt(tn)
power <- colMeans(tp)
power_se <- apply(tp, 2, sd)/sqrt(tn)
coverage_mean <- colMeans(coverage)
coverage_se <- apply(coverage, 2, sd)/sqrt(tn)
width_mean <- colMeans(width)
width_se <- apply(width, 2, sd)/sqrt(tn)

## plot

par(mfrow=c(2,2), cex.lab=1.5, cex.axis=1.3, mar=c(5,5,2,1), mgp=c(3,1,0))

plot(alpha_set, fwer, pch=16, ylim=c(0,0.22), xlab=expression(alpha), ylab="FWER")
arrows(alpha_set, fwer-fwer_se, alpha_set, fwer+fwer_se, length=0, angle=90, code=3, lwd=2)
segments(0.05, 0.05, 0.2, 0.2, lty=2, lwd=2)

plot(alpha_set, power, pch=16, ylim=c(0,1), xlab=expression(alpha), ylab="Power")
arrows(alpha_set, power-power_se, alpha_set, power+power_se, length=0, angle=90, code=3, lwd=2)

plot(1-alpha_set, coverage_mean, pch=16, ylim=c(0.7,1), xlab=expression(alpha), ylab="Rank interval coverage")
arrows(1-alpha_set, coverage_mean-coverage_se, 1-alpha_set, coverage_mean+coverage_se, length=0, angle=90, code=3, lwd=2)
segments(0.8, 0.8, 0.95, 0.95, lty=2, lwd=2)

plot(1-alpha_set, width_mean, pch=16, ylim=c(1,m), xlab=expression(alpha), ylab="Average rank interval width")
arrows(1-alpha_set, width_mean-width_se, 1-alpha_set, width_mean+width_se, length=0, angle=90, code=3, lwd=2)




## partitioning based method

p_X <- 20
sd_X <- 3
sd_Y <- 1
n_train <- 500

set.seed(12)
beta_1 <- runif(p_X)
beta_2 <- runif(p_X)
mu_X <- runif(p_X)

X_train <- mvrnorm(n_train, mu=mu_X, Sigma=sd_X*diag(p_X))
Y_train <- sapply(1:n_train, function(i) rnorm(1, mean=sum(beta_1*X_train[i,])+2*sin(sum(beta_2*X_train[i,])), sd=sd_Y*sqrt(sum(X_train[i,]^2))/p_X))
rf <- randomForest(X_train, Y_train)

n <- 1000
m <- 5
tn <- 500

alpha_set <- seq(0.05, 0.2, 0.025)
fd <- matrix(0, nrow=tn, ncol=length(alpha_set))
tp <- matrix(0, nrow=tn, ncol=length(alpha_set))
coverage <- matrix(0, nrow=tn, ncol=length(alpha_set))
width <- matrix(0, nrow=tn, ncol=length(alpha_set))

pb <- txtProgressBar(min=0, max=tn, style=3)

for(k in 1:tn) {
  X <- mvrnorm(n, mu=mu_X, Sigma=sd_X*diag(p_X))
  Y <- sapply(1:n, function(i) rnorm(1, mean=sum(beta_1*X[i,])+2*sin(sum(beta_2*X[i,])), sd=sd_Y*sqrt(sum(X[i,]^2))/p_X))
  muhat <- predict(rf, X)
  
  X_test <- mvrnorm(m, mu=mu_X, Sigma=sd_X*diag(p_X))
  Y_test <- sapply(1:m, function(i) rnorm(1, mean=sum(beta_1*X_test[i,])+2*sin(sum(beta_2*X_test[i,])), sd=sd_Y*sqrt(sum(X_test[i,]^2))/p_X))
  muhat_test <- predict(rf, X_test)
  
  test_pairs <- t(combn(m, 2))
  test_abs_diff <- apply(test_pairs, 1, function(a) abs(muhat_test[a[1]]-muhat_test[a[2]]))
  test_null <- apply(test_pairs, 1, function(a) (muhat_test[a[1]]-muhat_test[a[2]])*(Y_test[a[1]]-Y_test[a[2]]) <= 0)
  test_non_null <- 1-test_null
  
  I_hat <- rank(muhat_test, ties.method="random")
  true_rank <- rank(Y_test[order(muhat_test)])
  
  num_par <- floor(n/m)
  score <- rep(0, num_par)
  
  for(l in 1:num_par) {
    ind_l <- ((l-1)*m+1):(l*m)
    Y_l <- Y[ind_l]
    muhat_l <- muhat[ind_l]
    score[l] <- max(apply(test_pairs, 1, function(a) abs(muhat_l[a[1]]-muhat_l[a[2]])*((muhat_l[a[1]]-muhat_l[a[2]])*(Y_l[a[1]]-Y_l[a[2]]) <= 0)))
  }
  
  bound <- c(sort(score), Inf)[ceiling((1-alpha_set)*(num_par+1))]
  rej_pairs <- sapply(bound, function(t) test_abs_diff > t)
  fd[k,] <- apply(rej_pairs, 2, function(v) sum(v*test_null) > 0)
  tp[k,] <- apply(rej_pairs, 2, function(v) sum(v*test_non_null)/max(sum(test_non_null), 1))
  
  for(j in seq_along(alpha_set)) {
    rank_interval <- cbind(rep(1, m), rep(m, m))
    rej_pairs_ind <- test_pairs[rej_pairs[,j], , drop=FALSE]
    
    if(nrow(rej_pairs_ind) > 0) {
      rej_pairs_ind_ordered <- t(apply(rej_pairs_ind, 1, function(v) sort(I_hat[v])))
      for(t in 1:m) {
        lower_t <- sum(rej_pairs_ind_ordered[,2] == t)+1
        upper_t <- m-sum(rej_pairs_ind_ordered[,1] == t)
        rank_interval[t,] <- c(lower_t, upper_t)
      }
    }
    
    coverage[k,j] <- min((rank_interval[,1] <= true_rank)*(true_rank <= rank_interval[,2]))
    width[k,j] <- mean(rank_interval[,2]-rank_interval[,1]+1)
  }
  
  setTxtProgressBar(pb, k)
}
close(pb)

fwer <- colMeans(fd)
fwer_se <- apply(fd, 2, sd)/sqrt(tn)
power <- colMeans(tp)
power_se <- apply(tp, 2, sd)/sqrt(tn)

coverage_mean <- colMeans(coverage)
coverage_se <- apply(coverage, 2, sd)/sqrt(tn)
width_mean <- colMeans(width)
width_se <- apply(width, 2, sd)/sqrt(tn)

## plot
par(mfrow=c(2,2), cex.lab=1.5, cex.axis=1.3, mar=c(5,5,2,1), mgp=c(3,1,0))

plot(alpha_set, fwer, pch=16, ylim=c(0,0.22), xlab=expression(alpha), ylab="FWER")
arrows(alpha_set, fwer-fwer_se, alpha_set, fwer+fwer_se, length=0, angle=90, code=3, lwd=2)
segments(0.05, 0.05, 0.2, 0.2, lty=2, lwd=2)

plot(alpha_set, power, pch=16, ylim=c(0,1), xlab=expression(alpha), ylab="Power")
arrows(alpha_set, power-power_se, alpha_set, power+power_se, length=0, angle=90, code=3, lwd=2)

plot(1-alpha_set, coverage_mean, pch=16, ylim=c(0.7,1), xlab=expression(alpha), ylab="Rank interval coverage")
arrows(1-alpha_set, coverage_mean-coverage_se, 1-alpha_set, coverage_mean+coverage_se, length=0, angle=90, code=3, lwd=2)
segments(0.8, 0.8, 0.95, 0.95, lty=2, lwd=2)

plot(1-alpha_set, width_mean, pch=16, ylim=c(1,m), xlab=expression(alpha), ylab="Average rank interval width")
arrows(1-alpha_set, width_mean-width_se, 1-alpha_set, width_mean+width_se, length=0, angle=90, code=3, lwd=2)


