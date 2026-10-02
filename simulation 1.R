## functions
pBH <- function(p_values,alpha)
{
  ##inputs
  #p_values : vector of p values
  #alpha : target level
  
  ##output
  #rejection threshold
  
  p_sorted <- sort(p_values)
  M <- length(p_values)
  ind <- which(p_sorted <= (1:M)/M*alpha)
  if(length(ind)>0) p_sorted[max(ind)] else 0
}

eBH <- function(e_values,alpha)
{
  ##inputs
  #e_values : vector of e values
  #alpha : target level
  
  ##output
  #rejection threshold
  
  M <- length(e_values)
  e_sorted <- sort(e_values,decreasing=TRUE)
  ind <- which((1:M)*e_sorted/M >= 1/alpha)
  if(length(ind)>0) e_sorted[max(ind)] else Inf
}

##-------------------------------------------------------------------------------------------------

library(MASS)
library(randomForest)
library(pbapply)
library(parallel)

##parameters
p_X <- 20 #feature dimension
sd_X <- 3
sd_Y <- 1

##sample sizes
n_train <- 500
n <- 200
m <- 20

set.seed(1234)
beta_1 <- runif(p_X)
beta_2 <- runif(p_X)
mu_X <- runif(p_X)

##training data
X_train <- mvrnorm(n_train,mu=mu_X,Sigma=sd_X*diag(p_X))
mu_train <- as.vector(X_train%*%beta_1)+2*sin(as.vector(X_train%*%beta_2))
Y_train <- rnorm(n_train,mean=mu_train,sd=sd_Y*sqrt(rowSums(X_train^2))/p_X)   #setting 2
#Y_train <- sapply(1:n_train,function(i) rnorm(1,mean=sum(beta_1*X_train[i,]),sd=sd_Y))  #setting 1

##fitting
rf <- randomForest(X_train,Y_train,mtry=floor(p_X/2),ntree=1000)
#rlm <- lm(Y_train~I(X_train))

tn <- 500
alpha_set <- seq(0.05,0.2,0.025)
n_alpha <- length(alpha_set)
n_1 <- floor(n/2)

test_pairs <- t(combn(m,2))
n_pairs <- nrow(test_pairs)
pair_1 <- test_pairs[,1]
pair_2 <- test_pairs[,2]
pair_factor <- choose(m-2,2)/(n_1+1)
pairs_by_index <- lapply(1:m,function(a) which(pair_1==a | pair_2==a))

##main simulation
sim_one <- function(k)
{
  X <- mvrnorm(n,mu=mu_X,Sigma=sd_X*diag(p_X))
  mu_true <- as.vector(X%*%beta_1)+2*sin(as.vector(X%*%beta_2))
  Y <- rnorm(n,mean=mu_true,sd=sd_Y*sqrt(rowSums(X^2))/p_X)
  #Y <- sapply(1:n,function(i) rnorm(1,mean=sum(beta_1*X[i,]),sd=sd_Y))
  muhat <- predict(rf,X)
  #muhat <- predict(rlm,newdata=data.frame(X_train=I(X)))
  
  X_test <- mvrnorm(m,mu=mu_X,Sigma=sd_X*diag(p_X))
  mu_true_test <- as.vector(X_test%*%beta_1)+2*sin(as.vector(X_test%*%beta_2))
  Y_test <- rnorm(m,mean=mu_true_test,sd=sd_Y*sqrt(rowSums(X_test^2))/p_X)
  #Y_test <- sapply(1:m,function(i) rnorm(1,mean=sum(beta_1*X_test[i,]),sd=sd_Y))
  muhat_test <- predict(rf,newdata=X_test)
  #muhat_test <- predict(rlm,newdata=data.frame(X_train=I(X_test)))
  
  cal_diff <- muhat[1:n_1]-muhat[(n_1+1):(2*n_1)]
  cal_abs_diff <- abs(cal_diff)
  cal_Y_null <- cal_diff*(Y[1:n_1]-Y[(n_1+1):(2*n_1)])<=0
  cal_null_sorted <- sort(cal_abs_diff[cal_Y_null])
  n_cal_null <- length(cal_null_sorted)
  
  test_diff <- muhat_test[pair_1]-muhat_test[pair_2]
  test_abs_diff <- abs(test_diff)
  test_null <- test_diff*(Y_test[pair_1]-Y_test[pair_2])<=0
  test_non_null <- !test_null
  pos <- sum(test_non_null)
  true_rank <- rank(Y_test[order(muhat_test)])
  
  t_set <- sort(test_abs_diff)
  total_count <- n_pairs-findInterval(t_set,t_set)
  numer_set <- n_cal_null-findInterval(t_set,cal_null_sorted)+1
  
  index_count <- matrix(0L,nrow=m,ncol=n_pairs)
  for(a in 1:m)
  {
    vals <- sort(test_abs_diff[pairs_by_index[[a]]])
    index_count[a,] <- length(vals)-findInterval(t_set,vals)
  }
  
  T_rej_mat <- matrix(Inf,nrow=n_pairs,ncol=n_alpha)
  for(l in seq_len(n_pairs))
  {
    denom <- pmax(total_count-index_count[pair_1[l],]-index_count[pair_2[l],]+as.integer(test_abs_diff[l]>t_set),1)
    FDP_l <- numer_set/denom*pair_factor
    for(j in seq_len(n_alpha))
    {
      ind <- which(FDP_l<=alpha_set[j])
      if(length(ind)>0) T_rej_mat[l,j] <- t_set[ind[1]]
    }
  }
  
  p_val <- (n_cal_null-findInterval(test_abs_diff,cal_null_sorted)+1)/(n_1+1)
  out <- matrix(0,nrow=n_alpha,ncol=18)
  
  get_rank_result <- function(rej_pairs)
  {
    rank_interval <- cbind(rep(1,m),rep(m,m))
    if(sum(rej_pairs)>0)
    {
      I_hat <- rank(muhat_test)
      rej_idx <- which(rej_pairs)
      r1 <- I_hat[pair_1[rej_idx]]
      r2 <- I_hat[pair_2[rej_idx]]
      low <- pmin(r1,r2)
      high <- pmax(r1,r2)
      rank_interval[,1] <- tabulate(high,nbins=m)+1
      rank_interval[,2] <- m-tabulate(low,nbins=m)
    }
    risk_upper <- sum(pmax(true_rank-rank_interval[,2],0))/max(sum(m-rank_interval[,2]),1)
    risk_lower <- sum(pmax(rank_interval[,1]-true_rank,0))/max(sum(rank_interval[,1]-1),1)
    coverage <- mean(rank_interval[,1]<=true_rank & true_rank<=rank_interval[,2])
    width <- mean(rank_interval[,2]-rank_interval[,1]+1)
    c(risk_upper,risk_lower,coverage,width)
  }
  
  for(j in seq_len(n_alpha))
  {
    alpha <- alpha_set[j]
    T_rej <- T_rej_mat[,j]
    finite_T <- is.finite(T_rej)
    e_denom <- rep(1,n_pairs)
    e_denom[finite_T] <- n_cal_null-findInterval(T_rej[finite_T],cal_null_sorted)+1
    e_val <- as.numeric(test_abs_diff>T_rej)*(n_1+1)/e_denom
    U <- runif(1)
    e_val_U <- e_val/U
    
    rej_pairs <- e_val>=eBH(e_val,alpha)
    rej_pairs_U <- e_val_U>=eBH(e_val_U,alpha)
    rej_pairs_p <- p_val<=pBH(p_val,alpha)
    
    R <- sum(rej_pairs)
    V <- sum(rej_pairs & test_null)
    tr_pos <- sum(rej_pairs & test_non_null)
    R_U <- sum(rej_pairs_U)
    V_U <- sum(rej_pairs_U & test_null)
    tr_pos_U <- sum(rej_pairs_U & test_non_null)
    R_p <- sum(rej_pairs_p)
    V_p <- sum(rej_pairs_p & test_null)
    tr_pos_p <- sum(rej_pairs_p & test_non_null)
    
    rank_res <- get_rank_result(rej_pairs)
    rank_res_U <- get_rank_result(rej_pairs_U)
    rank_res_p <- get_rank_result(rej_pairs_p)
    
    out[j,] <- c(V/max(1,R),tr_pos/max(1,pos),V_U/max(1,R_U),tr_pos_U/max(1,pos),V_p/max(1,R_p),tr_pos_p/max(1,pos),rank_res,rank_res_U,rank_res_p)
  }
  out
}

cl <- makeCluster(max(1,detectCores()-1))
clusterEvalQ(cl,{library(MASS);library(randomForest)})
clusterExport(cl,c("pBH","eBH","p_X","sd_X","sd_Y","n","m","n_1","mu_X","beta_1","beta_2","alpha_set","n_alpha","n_pairs","pair_1","pair_2","pair_factor","pairs_by_index","rlm","rf"),envir=environment())

res_list <- pblapply(seq_len(tn),sim_one,cl=cl)
stopCluster(cl)

res_array <- simplify2array(res_list)

fdp <- t(res_array[,1,])
tpp <- t(res_array[,2,])
fdp_U <- t(res_array[,3,])
tpp_U <- t(res_array[,4,])
fdp_p <- t(res_array[,5,])
tpp_p <- t(res_array[,6,])

risk_upper <- t(res_array[,7,])
risk_lower <- t(res_array[,8,])
coverage <- t(res_array[,9,])
width <- t(res_array[,10,])

risk_upper_U <- t(res_array[,11,])
risk_lower_U <- t(res_array[,12,])
coverage_U <- t(res_array[,13,])
width_U <- t(res_array[,14,])

risk_upper_p <- t(res_array[,15,])
risk_lower_p <- t(res_array[,16,])
coverage_p <- t(res_array[,17,])
width_p <- t(res_array[,18,])

fdr <- colMeans(fdp)
fdr_se <- apply(fdp,2,sd)/sqrt(tn)
power <- colMeans(tpp)
power_se <- apply(tpp,2,sd)/sqrt(tn)

fdr_U <- colMeans(fdp_U)
fdr_U_se <- apply(fdp_U,2,sd)/sqrt(tn)
power_U <- colMeans(tpp_U)
power_U_se <- apply(tpp_U,2,sd)/sqrt(tn)

fdr_p <- colMeans(fdp_p)
fdr_p_se <- apply(fdp_p,2,sd)/sqrt(tn)
power_p <- colMeans(tpp_p)
power_p_se <- apply(tpp_p,2,sd)/sqrt(tn)

risk_upper_mean <- colMeans(risk_upper)
risk_upper_se <- apply(risk_upper,2,sd)/sqrt(tn)
risk_lower_mean <- colMeans(risk_lower)
risk_lower_se <- apply(risk_lower,2,sd)/sqrt(tn)
coverage_mean <- colMeans(coverage)
coverage_se <- apply(coverage,2,sd)/sqrt(tn)
width_mean <- colMeans(width)
width_se <- apply(width,2,sd)/sqrt(tn)

risk_upper_U_mean <- colMeans(risk_upper_U)
risk_upper_U_se <- apply(risk_upper_U,2,sd)/sqrt(tn)
risk_lower_U_mean <- colMeans(risk_lower_U)
risk_lower_U_se <- apply(risk_lower_U,2,sd)/sqrt(tn)
coverage_U_mean <- colMeans(coverage_U)
coverage_U_se <- apply(coverage_U,2,sd)/sqrt(tn)
width_U_mean <- colMeans(width_U)
width_U_se <- apply(width_U,2,sd)/sqrt(tn)

risk_upper_p_mean <- colMeans(risk_upper_p)
risk_upper_p_se <- apply(risk_upper_p,2,sd)/sqrt(tn)
risk_lower_p_mean <- colMeans(risk_lower_p)
risk_lower_p_se <- apply(risk_lower_p,2,sd)/sqrt(tn)
coverage_p_mean <- colMeans(coverage_p)
coverage_p_se <- apply(coverage_p,2,sd)/sqrt(tn)
width_p_mean <- colMeans(width_p)
width_p_se <- apply(width_p,2,sd)/sqrt(tn)

res <- cbind(fdr,fdr_se,fdr_U,fdr_U_se,fdr_p,fdr_p_se,
             power,power_se,power_U,power_U_se,power_p,power_p_se,
             risk_upper_mean,risk_upper_se,risk_upper_U_mean,risk_upper_U_se,risk_upper_p_mean,risk_upper_p_se,
             risk_lower_mean,risk_lower_se,risk_lower_U_mean,risk_lower_U_se,risk_lower_p_mean,risk_lower_p_se,
             coverage_mean,coverage_se,coverage_U_mean,coverage_U_se,coverage_p_mean,coverage_p_se,
             width_mean,width_se,width_U_mean,width_U_se,width_p_mean,width_p_se)

#plot
col1 <- "#0072B2"
col2 <- "#D55E00"
col3 <- "#CC79A7"
pch_set <- c(16,17,4)

par(mfrow=c(3,2),cex.lab=2,cex.axis=1.7,mar=c(5,6,2,0),xpd=NA,oma=c(0,0,0,4),mgp=c(3,1,0))

matplot(alpha_set,res[,c(1,3,5)],pch=pch_set,type="b",ylim=c(0,0.22),xlab=expression(alpha),ylab="FDR",cex=0.8,col=c(col1,col3,col2))
arrows(alpha_set,res[,1]-res[,2],alpha_set,res[,1]+res[,2],length=0,angle=90,code=3,lwd=2,col=col1)
arrows(alpha_set,res[,3]-res[,4],alpha_set,res[,3]+res[,4],length=0,angle=90,code=3,lwd=2,col=col3)
arrows(alpha_set,res[,5]-res[,6],alpha_set,res[,5]+res[,6],length=0,angle=90,code=3,lwd=2,col=col2)
segments(0.05,0.05,0.2,0.2,lty=2,lwd=2)
legend("topleft",legend=c("e-BH","U-eBH","p-BH"),pch=pch_set,col=c(col1,col3,col2),bty="n",cex=1.2)

matplot(alpha_set,res[,c(7,9,11)],pch=pch_set,type="b",ylim=c(0,1),xlab=expression(alpha),ylab="Power",cex=0.8,col=c(col1,col3,col2))
arrows(alpha_set,res[,7]-res[,8],alpha_set,res[,7]+res[,8],length=0,angle=90,code=3,lwd=2,col=col1)
arrows(alpha_set,res[,9]-res[,10],alpha_set,res[,9]+res[,10],length=0,angle=90,code=3,lwd=2,col=col3)
arrows(alpha_set,res[,11]-res[,12],alpha_set,res[,11]+res[,12],length=0,angle=90,code=3,lwd=2,col=col2)

matplot(alpha_set,res[,c(13,15,17)],pch=pch_set,type="b",ylim=c(0,0.22),xlab=expression(alpha),ylab="Risk+",cex=0.8,col=c(col1,col3,col2))
arrows(alpha_set,res[,13]-res[,14],alpha_set,res[,13]+res[,14],length=0,angle=90,code=3,lwd=2,col=col1)
arrows(alpha_set,res[,15]-res[,16],alpha_set,res[,15]+res[,16],length=0,angle=90,code=3,lwd=2,col=col3)
arrows(alpha_set,res[,17]-res[,18],alpha_set,res[,17]+res[,18],length=0,angle=90,code=3,lwd=2,col=col2)
segments(0.05,0.05,0.2,0.2,lty=2,lwd=2)

matplot(alpha_set,res[,c(19,21,23)],pch=pch_set,type="b",ylim=c(0,0.22),xlab=expression(alpha),ylab="Risk-",cex=0.8,col=c(col1,col3,col2))
arrows(alpha_set,res[,19]-res[,20],alpha_set,res[,19]+res[,20],length=0,angle=90,code=3,lwd=2,col=col1)
arrows(alpha_set,res[,21]-res[,22],alpha_set,res[,21]+res[,22],length=0,angle=90,code=3,lwd=2,col=col3)
arrows(alpha_set,res[,23]-res[,24],alpha_set,res[,23]+res[,24],length=0,angle=90,code=3,lwd=2,col=col2)
segments(0.05,0.05,0.2,0.2,lty=2,lwd=2)

matplot(1-alpha_set,res[,c(25,27,29)],pch=pch_set,type="b",ylim=c(0,1),xlab=expression(1-alpha),ylab="Mean coverage rate",cex=0.8,col=c(col1,col3,col2))
arrows(1-alpha_set,res[,25]-res[,26],1-alpha_set,res[,25]+res[,26],length=0,angle=90,code=3,lwd=2,col=col1)
arrows(1-alpha_set,res[,27]-res[,28],1-alpha_set,res[,27]+res[,28],length=0,angle=90,code=3,lwd=2,col=col3)
arrows(1-alpha_set,res[,29]-res[,30],1-alpha_set,res[,29]+res[,30],length=0,angle=90,code=3,lwd=2,col=col2)
segments(0.8,0.8,0.95,0.95,lty=2,lwd=2)

matplot(alpha_set,res[,c(31,33,35)],pch=pch_set,type="b",ylim=c(0,m),xlab=expression(alpha),ylab="Mean width",cex=0.8,col=c(col1,col3,col2))
arrows(alpha_set,res[,31]-res[,32],alpha_set,res[,31]+res[,32],length=0,angle=90,code=3,lwd=2,col=col1)
arrows(alpha_set,res[,33]-res[,34],alpha_set,res[,33]+res[,34],length=0,angle=90,code=3,lwd=2,col=col3)
arrows(alpha_set,res[,35]-res[,36],alpha_set,res[,35]+res[,36],length=0,angle=90,code=3,lwd=2,col=col2)


