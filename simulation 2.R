library(MASS)
library(randomForest)

#parameters
p_X <- 20

sd_X <- 3
sd_Y <- 1

n_train <- 500

set.seed(12)
beta_1 <- runif(p_X)
beta_2 <- runif(p_X)
mu_X <- runif(p_X)

#generate data
X_train <- mvrnorm(n_train, mu=mu_X, Sigma=sd_X*diag(p_X))
Y_train <- sapply(1:n_train, function(i) rnorm(1,mean=sum(beta_1*X_train[i,])+2*sin(sum(beta_2*X_train[i,])), sd=sd_Y*sqrt(sum(X_train[i,]^2))/p_X))


rf <- randomForest(X_train,Y_train)

n <- 1000
m <- 5

tn <- 500

alpha_set <- seq(0.05,0.2,0.025)
fd <- matrix(0,nrow=tn,ncol=length(alpha_set))
tp <- matrix(0,nrow=tn,ncol=length(alpha_set))

pb <- txtProgressBar(min = 0, max = tn, style = 3)
for(k in 1:tn)
{
  
  X <- mvrnorm(n, mu=mu_X, Sigma=sd_X*diag(p_X))
  Y <- sapply(1:n, function(i) rnorm(1,mean=sum(beta_1*X[i,])+2*sin(sum(beta_2*X[i,])), sd=sd_Y*sqrt(sum(X[i,]^2))/p_X))
  
  muhat <- predict(rf,X)
  
  X_test <- mvrnorm(m, mu=mu_X, Sigma=sd_X*diag(p_X))
  Y_test <- sapply(1:m, function(i) rnorm(1,mean=sum(beta_1*X_test[i,])+2*sin(sum(beta_2*X_test[i,])), sd=sd_Y*sqrt(sum(X_test[i,]^2))/p_X))
  
  muhat_test <- predict(rf,X_test)
  
  n_1 <- floor(n/2)
  cal_abs_diff <- abs(muhat[1:n_1]-muhat[(n_1+1):(2*n_1)])
  cal_Y_null <- ( (muhat[1:n_1]-muhat[(n_1+1):(2*n_1)])*(Y[1:n_1]-Y[(n_1+1):(2*n_1)]) <= 0 )
 
  test_pairs <- t(combn(m,2))
  test_abs_diff <- apply(test_pairs,1, function(a) abs(muhat_test[a[1]] - muhat_test[a[2]]))
  test_null <- apply(test_pairs,1, function(a) (muhat_test[a[1]] - muhat_test[a[2]])*(Y_test[a[1]]-Y_test[a[2]]) <= 0)
  test_non_null <- 1-test_null
  
  
  p_val <- sapply(test_abs_diff, function(v) (sum((v < cal_abs_diff)*cal_Y_null)+1)/(n_1+1))
  
  for(j in 1:length(alpha_set))
  {
    alpha <- alpha_set[j]
    rej_pairs <- (p_val <= (2*alpha/(m*(m-1))))
    fd[k,j] <- (sum(rej_pairs*test_null) > 0)
    tp[k,j] <- sum(rej_pairs*test_non_null)/max(sum(test_non_null),1)
  
  }
  setTxtProgressBar(pb, k)
}

fwer <- colMeans(fd)
fwer_se <- apply(fd,2,sd)/sqrt(tn)
power <- colMeans(tp)
power_se <- apply(tp,2,sd)/sqrt(tn)


#plot

par(mfrow=c(1,2), cex.lab = 2, cex.axis=1.7,mar=c(5,6,2,0), xpd=NA, oma = c(0, 0, 0, 4), mgp = c(3, 1, 0))
plot(alpha_set,fwer,pch=16,ylim=c(0,0.22),xlab=expression(alpha),ylab="FWER")
arrows(alpha_set, fwer-fwer_se, alpha_set, fwer+fwer_se, length=0, angle=90, code=3,lwd=2)
segments(0.05,0.05,0.2,0.2,lty=2,lwd=2)
plot(alpha_set,power,pch=16,ylim=c(0,0.5),xlab=expression(alpha),ylab="Power")
arrows(alpha_set, power-power_se, alpha_set, power+power_se, length=0, angle=90, code=3,lwd=2)




