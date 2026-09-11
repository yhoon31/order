##functions

#BH procedure
pBH <- function(p_values, alpha) 
{
  ##inputs
  #p_values : vector of p values
  #alpha : target level
  
  ##output
  #rejection threshold
  
  p_sorted <- sort(p_values)
  m <- length(p_values)
  if(min(p_sorted*m/(1:m)) <= alpha)
  {
    threshold <- p_sorted[max(which(p_sorted <= (1:m) / m * alpha))]
  }
  else
  {
    threshold <- 0
  }
  return(threshold)
  
}

#e-BH procedure
eBH <- function(e_values,alpha)
{
  ##inputs
  #e_values : vector of e values
  #alpha : target level
  
  ##output
  #rejection threshold
  
  m <- length(e_values)
  e_sorted <- sort(e_values,decreasing=TRUE)
  e_adjusted <- (1:m)*e_sorted/m
  if(max(e_adjusted) < 1/alpha)
  {
    threshold <- Inf
  }  else
  {
    threshold <- e_sorted[max(which(e_adjusted >= 1/alpha))]
  }
  return(threshold)
}


##-------------------------------------------------------------------------------------------------

library(MASS)
library(randomForest)
library(glmnet)

##parameters
p_X <- 20 #feature dimension
sd_X <- 3 
sd_Y <- 1

set.seed(1234)
beta_1 <- runif(p_X)
beta_2 <- runif(p_X)
mu_X <- runif(p_X)

##sample sizes
n_train <- 500
n <- 200
m <- 20

##training data
X_train <- mvrnorm(n_train, mu=mu_X, Sigma=sd_X*diag(p_X))
Y_train <- sapply(1:n_train, function(i) rnorm(1,mean=sum(beta_1*X_train[i,])+2*sin(sum(beta_2*X_train[i,])), sd=sd_Y*sqrt(sum(X_train[i,]^2))/p_X)) #Setting 2
#Y_train <- sapply(1:n_train, function(i) rnorm(1,mean=sum(beta_1*X_train[i,]), sd=sd_Y)) #Setting 1

##fitting
rf <- randomForest(X_train,Y_train,mtry=floor(ncol(X_train)/2),ntree=1000) #random forest
#rlm <- lm(Y_train~I(X_train)) #linear regression

tn <- 500

alpha_set <- seq(0.05,0.2,0.025)
fdp <- matrix(0,nrow=tn,ncol=length(alpha_set))
tpp <- matrix(0,nrow=tn,ncol=length(alpha_set))
fdp_U <- matrix(0,nrow=tn,ncol=length(alpha_set))
tpp_U <- matrix(0,nrow=tn,ncol=length(alpha_set))
fdp_p <- matrix(0,nrow=tn,ncol=length(alpha_set))
tpp_p <- matrix(0,nrow=tn,ncol=length(alpha_set))
risk_upper <- matrix(0,nrow=tn,ncol=length(alpha_set))
risk_lower <- matrix(0,nrow=tn,ncol=length(alpha_set))
coverage <- matrix(0,nrow=tn,ncol=length(alpha_set))
width <- matrix(0,nrow=tn,ncol=length(alpha_set))

##main simulation
pb <- txtProgressBar(min = 0, max = tn, style = 3)
for(k in 1:tn)
{
  ##data generation
  X <- mvrnorm(n, mu=mu_X, Sigma=sd_X*diag(p_X))
  Y <- sapply(1:n, function(i) rnorm(1,mean=sum(beta_1*X[i,])+2*sin(sum(beta_2*X[i,])), sd=sd_Y*sqrt(sum(X[i,]^2))/p_X))
  #Y <- sapply(1:n, function(i) rnorm(1,mean=sum(beta_1*X[i,]), sd=sd_Y))
  
  muhat <- predict(rf,X)
  #muhat <- predict(rlm, newdata = data.frame(X_train = I(X)))
 
  X_test <- mvrnorm(m, mu=mu_X, Sigma=sd_X*diag(p_X))
  Y_test <- sapply(1:m, function(i) rnorm(1,mean=sum(beta_1*X_test[i,])+2*sin(sum(beta_2*X_test[i,])), sd=sd_Y*sqrt(sum(X_test[i,]^2))/p_X))
  #Y_test <- sapply(1:m, function(i) rnorm(1,mean=sum(beta_1*X_test[i,]), sd=sd_Y))
  
  muhat_test <- predict(rf,newdata=X_test)
  #muhat_test <- predict(rlm, newdata = data.frame(X_train = I(X_test)))

  
  ##compute statistics
  n_1 <- floor(n/2)
  cal_abs_diff <- abs(muhat[1:n_1]-muhat[(n_1+1):(2*n_1)])
  cal_Y_null <- ( (muhat[1:n_1]-muhat[(n_1+1):(2*n_1)])*(Y[1:n_1]-Y[(n_1+1):(2*n_1)]) <= 0 )
  
  test_pairs <- t(combn(m,2))
  test_abs_diff <- apply(test_pairs,1, function(a) abs(muhat_test[a[1]] - muhat_test[a[2]]))
  test_null <- apply(test_pairs,1, function(a) (muhat_test[a[1]] - muhat_test[a[2]])*(Y_test[a[1]]-Y_test[a[2]]) <= 0)
  test_non_null <- 1-test_null
  
  true_rank <- rank(Y_test[order(muhat_test)])
  
  FDP_hat <- function(t,j)
  {
    numer <- sum((cal_abs_diff > t)*cal_Y_null)+1
    denom <- max(1,sum(test_abs_diff > t))
    FDP_t <- (numer/denom)*(m*(m-1))/(2*(n_1+1))
    return(FDP_t)
  }
  FDP_hat <- Vectorize(FDP_hat)
  
  t_set <- sort(test_abs_diff)
  FDP_set <- FDP_hat(t_set)
  
  ##run procedures at different levels
  for(j in 1:length(alpha_set))
  {
    alpha <- alpha_set[j]
    alpha_tilde <- alpha
    
    if(min(FDP_set) <= alpha_tilde)
    {
      T_rej <- t_set[min(which(FDP_set <= alpha_tilde))]
    } else{
      T_rej <- Inf
    } 
    
    U <- runif(1)
    
    e_val <- ((test_abs_diff > T_rej)/(sum((cal_abs_diff>T_rej)*cal_Y_null)+1))*(n_1+1)
    p_val <- sapply(test_abs_diff, function(v) (sum((v > cal_abs_diff)*cal_Y_null)+1)/(n_1+1))
    
    e_val_U <- ((test_abs_diff > T_rej)/(sum((cal_abs_diff>T_rej)*cal_Y_null)+1))*(n_1+1)/U
    
    rej_pairs <- (e_val >= eBH(e_val,alpha))
    rej_pairs_U <- (e_val_U >= eBH(e_val,alpha))
    
    #eBH FDP/power
    R <- sum(rej_pairs)
    V <- sum(rej_pairs*test_null)
    pos <- sum(test_non_null)
    tr_pos <- sum(rej_pairs*test_non_null)
    
    fdp[k,j] <- V/max(1,R)
    tpp[k,j] <- tr_pos/max(1,pos)
    
    #U-eBH FDP/power
    R_U <- sum(rej_pairs_U)
    V_U <- sum(rej_pairs_U*test_null)
    tr_pos_U <- sum(rej_pairs_U*test_non_null)
    
    fdp_U[k,j] <- V_U/max(1,R_U)
    tpp_U[k,j] <- tr_pos_U/max(1,pos)
    
    #pBH FDP/power
    
    rej_pairs_p <- (p_val <= pBH(p_val,alpha))
    
    R_p <- sum(rej_pairs_p)
    V_p <- sum(rej_pairs_p*test_null)
    tr_pos_p <- sum(rej_pairs_p*test_non_null)
    
    fdp_p[k,j] <- V_p/max(1,R_p)
    tpp_p[k,j] <- tr_pos_p/max(1,pos)
    
    
    #rank interval
    rank_interval <- cbind(rep(1,m),rep(m,m))
    
    if(sum(rej_pairs)!=0)
    {
      rej_pairs_ind <-  matrix(test_pairs[rej_pairs==T,],ncol=2)
      
      I_hat <- rank(muhat_test)
      
      rej_pairs_ind_ordered <- t(apply(rej_pairs_ind,1, function(v) sort(c(I_hat[v[1]],I_hat[v[2]])) ))
      
      rej_pairs_ind_ordered <- rej_pairs_ind_ordered[order(rej_pairs_ind_ordered[,1],rej_pairs_ind_ordered[,2]),]
      
      for(t in 1:m)
      {
        lower_t <- sum(rej_pairs_ind_ordered[,2]==t)+1
        upper_t <- m-sum(rej_pairs_ind_ordered[,1]==t)
        rank_interval[t,] <- c(lower_t,upper_t)
      }
    }
    
    risk_upper[k,j] <- sum(pmax(true_rank - rank_interval[,2],0))/max(sum(m-rank_interval[,2]),1)
    risk_lower[k,j] <- sum(pmax(rank_interval[,1]-true_rank,0))/max(sum(rank_interval[,1]-1),1)
    coverage[k,j] <- mean((rank_interval[,1] <= true_rank)*(true_rank <= rank_interval[,2]))
    width[k,j] <- mean(rank_interval[,2]-rank_interval[,1]+1)
    
  }
  
  
  setTxtProgressBar(pb, k)
}

##compute fdr/power/risk
fdr <- colMeans(fdp)
fdr_se <- apply(fdp,2,sd)/sqrt(tn)
fdr_U <- colMeans(fdp_U)
fdr_U_se <- apply(fdp_U,2,sd)/sqrt(tn)
fdr_p <- colMeans(fdp_p)
fdr_p_se <- apply(fdp_p,2,sd)/sqrt(tn)
power <- colMeans(tpp)
power_se <- apply(tpp,2,sd)/sqrt(tn)
power_p <- colMeans(tpp_p)
power_p_se <- apply(tpp_p,2,sd)/sqrt(tn)
power_U <- colMeans(tpp_U)
power_U_se <- apply(tpp_U,2,sd)/sqrt(tn)

risk_upper_mean <- colMeans(risk_upper)
risk_upper_se <- apply(risk_upper,2,sd)/sqrt(tn)
risk_lower_mean <- colMeans(risk_lower)
risk_lower_se <- apply(risk_lower,2,sd)/sqrt(tn)

coverage_mean <- colMeans(coverage)
coverage_se <- apply(coverage,2,sd)/sqrt(tn)
width_mean <- colMeans(width)
width_se <- apply(width,2,sd)/sqrt(tn)

res <- cbind(fdr,fdr_se,power,power_se,risk_upper_mean,risk_upper_se,risk_lower_mean,risk_lower_se,coverage_mean,coverage_se,width_mean,width_se,fdr_U,fdr_U_se,power_U,power_U_se,fdr_p,fdr_p_se,power_p,power_p_se)

##plot
par(mfrow=c(3,2), cex.lab = 2, cex.axis=1.7,mar=c(5,6,2,0), xpd=NA, oma = c(0, 0, 0, 4), mgp = c(3, 1, 0))
plot(alpha_set,fdr,pch=16,ylim=c(0,0.22),xlab=expression(alpha),ylab="FDR",cex=0.8)
arrows(alpha_set, fdr-fdr_se, alpha_set, fdr+fdr_se, length=0, angle=90, code=3,lwd=2)
segments(0.05,0.05,0.2,0.2,lty=2,lwd=2)
plot(alpha_set,power,pch=16,ylim=c(0,1),xlab=expression(alpha),ylab="Power",cex=0.8)
arrows(alpha_set, power-power_se, alpha_set, power+power_se, length=0, angle=90, code=3,lwd=2)

plot(alpha_set,risk_upper_mean ,pch=16,ylim=c(0,0.22),xlab=expression(alpha),ylab="Risk+",cex=0.8)
arrows(alpha_set, risk_upper_mean - risk_upper_se, alpha_set, risk_upper_mean + risk_upper_se, length=0, angle=90, code=3,lwd=2)
segments(0.05,0.05,0.2,0.2,lty=2,lwd=2)
plot(alpha_set,risk_lower_mean ,pch=16,ylim=c(0,0.22),xlab=expression(alpha),ylab="Risk-",cex=0.8)
arrows(alpha_set, risk_lower_mean - risk_lower_se, alpha_set, risk_lower_mean + risk_lower_se, length=0, angle=90, code=3,lwd=2)
segments(0.05,0.05,0.2,0.2,lty=2,lwd=2)

plot(1-alpha_set,coverage_mean,pch=16,ylim=c(0,1),xlab=expression(1-alpha),ylab="Mean coverage rate",cex=0.8)
arrows(1-alpha_set, coverage_mean-coverage_se, 1-alpha_set, coverage_mean+coverage_se, length=0, angle=90, code=3,lwd=2)
plot(1-alpha_set,width_mean,pch=16,ylim=c(0,m),xlab=expression(alpha),ylab="Mean width",cex=0.8)
arrows(1-alpha_set, width_mean-width_se, 1-alpha_set, width_mean+width_se, length=0, angle=90, code=3,lwd=2)

