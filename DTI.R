##e-BH procedure
eBH <- function(e_values,alpha)
{
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

##----------------------------------------------------------------------

#load data
y_cal <- read.table("y_cal.txt")$V1
y_test_all <- read.table("y_test.txt")$V1
muhat <- read.table("muhat_cal.txt")$V1
muhat_test_all <- read.table("muhat_test.txt")$V1

set.seed(1)
perm <- sample(1:length(y_test_all)) #shuffle test set
y_test_all <- y_test_all[perm]
muhat_test_all <- muhat_test_all[perm]


##FWER control (Appendix)

n <- 10000
m <- 5
num_test_set <- 500

alpha_set <- seq(0.05,0.2,0.025)
fd <- matrix(0,nrow=num_test_set,ncol=length(alpha_set))
tp <- matrix(0,nrow=num_test_set,ncol=length(alpha_set))


n_1 <- floor(n/2)
cal_abs_diff <- abs(muhat[1:n_1]-muhat[(n_1+1):(2*n_1)])
cal_Y_null <- ( (muhat[1:n_1]-muhat[(n_1+1):(2*n_1)])*(y_cal[1:n_1]-y_cal[(n_1+1):(2*n_1)]) <= 0 )
test_pairs <- t(combn(m,2))


pb <- txtProgressBar(min = 0, max = num_test_set, style = 3)

for(t in 1:num_test_set)
{
  muhat_test <- abs(muhat_test_all[(m*(t-1)+1):(m*t)])
  Y_test <- y_test_all[(m*(t-1)+1):(m*t)]
  test_abs_diff <- apply(test_pairs,1, function(a) abs(muhat_test[a[1]] - muhat_test[a[2]]))
  
  test_null <- apply(test_pairs,1, function(a) (muhat_test[a[1]] - muhat_test[a[2]])*(Y_test[a[1]]-Y_test[a[2]]) <= 0)
  test_non_null <- 1-test_null
  
  
  p_val <- sapply(test_abs_diff, function(v) (sum((v < cal_abs_diff)*cal_Y_null)+1)/(n_1+1))
  
  for(j in 1:length(alpha_set))
  {
    alpha <- alpha_set[j]
    rej_pairs <- (p_val <= (2*alpha/(m*(m-1))))
    fd[t,j] <- (sum(rej_pairs*test_null) > 0)
    tp[t,j] <- sum(rej_pairs*test_non_null)/max(sum(test_non_null),1)
    
  }
  
  setTxtProgressBar(pb, t)
}


fwer <- colMeans(fd)
fwer_se <- apply(fd,2,sd)/sqrt(num_test_set)
power <- colMeans(tp)
power_se <- apply(tp,2,sd)/sqrt(num_test_set)


#plot
par(mfrow=c(1,2), cex.lab = 2, cex.axis=1.7,mar=c(5,6,2,0), xpd=NA, oma = c(0, 0, 0, 4), mgp = c(3, 1, 0))
plot(alpha_set,fwer,pch=16,ylim=c(0,0.22),xlab=expression(alpha),ylab="FWER")
arrows(alpha_set, fwer-fwer_se, alpha_set, fwer+fwer_se, length=0, angle=90, code=3,lwd=2)
segments(0.05,0.05,0.2,0.2,lty=2,lwd=2)
plot(alpha_set,power,pch=16,ylim=c(0,0.4),xlab=expression(alpha),ylab="Power")
arrows(alpha_set, power-power_se, alpha_set, power+power_se, length=0, angle=90, code=3,lwd=2)

#Figure

setwd("C:/Wharton/Inference on the order")
res1 <- read.table("DTI_FWER_5.txt")
res2 <- read.table("DTI_FWER_10.txt")
res3 <- read.table("DTI_FWER_20.txt")
alpha_set <- seq(0.05,0.2,0.025)

par(mfrow=c(2,3), cex.lab = 2.5, cex.axis=2,cex.main=2.5,font.main=1,mar=c(5,6,1,0), xpd=NA, oma = c(0, 0, 2.5, 4), mgp = c(3, 1, 0))
plot(alpha_set,res1[,1],pch=16,ylim=c(0,0.22),xlab="",ylab="FWER",cex=0.8)
title(main="Test size = 5", line=1.25)
arrows(alpha_set, res1[,1]-res1[,2], alpha_set, res1[,1]+res1[,2], length=0, angle=90, code=3,lwd=1.5)
segments(0.05,0.05,0.2,0.2,lty=2,lwd=2)
plot(alpha_set,res2[,1],pch=16,ylim=c(0,0.22),xlab="",ylab="",cex=0.8)
title(main="Test size = 10", line=1.25)
arrows(alpha_set, res2[,1]-res2[,2], alpha_set, res2[,1]+res2[,2], length=0, angle=90, code=3,lwd=1.5)
segments(0.05,0.05,0.2,0.2,lty=2,lwd=2)
plot(alpha_set,res3[,1],pch=16,ylim=c(0,0.22),xlab="",ylab="",cex=0.8)
title(main="Test size = 20", line=1.25)
arrows(alpha_set, res3[,1]-res3[,2], alpha_set, res3[,1]+res3[,2], length=0, angle=90, code=3,lwd=1.5)
segments(0.05,0.05,0.2,0.2,lty=2,lwd=2)

plot(alpha_set,res1[,3],pch=16,ylim=c(0,0.4),xlab=expression(alpha),ylab="Power",cex=0.8)
arrows(alpha_set, res1[,3]-res1[,4], alpha_set, res1[,3]+res1[,4], length=0, angle=90, code=3,lwd=1.5)
plot(alpha_set,res2[,3],pch=16,ylim=c(0,0.4),xlab=expression(alpha),ylab="",cex=0.8)
arrows(alpha_set, res2[,3]-res2[,4], alpha_set, res2[,3]+res2[,4], length=0, angle=90, code=3,lwd=1.5)
plot(alpha_set,res3[,3],pch=16,ylim=c(0,0.4),xlab=expression(alpha),ylab="",cex=0.8)
arrows(alpha_set, res3[,3]-res3[,4], alpha_set, res3[,3]+res3[,4], length=0, angle=90, code=3,lwd=1.5)

#1080x610



##FDR control (main paper)

y_cal_all <- read.table("y_cal.txt")$V1
y_test_all <- read.table("y_test.txt")$V1
muhat_cal_all <- read.table("muhat_cal.txt")$V1
muhat_test_all <- read.table("muhat_test.txt")$V1

n <- 200
m <- 20
num_test_set <- 80

set.seed(1234)
perm_cal <- sample(1:length(y_cal_all)) #shuffle calibration set
y_cal_all <- y_cal_all[perm_cal]
muhat_cal_all <- muhat_cal_all[perm_cal]
perm <- sample(1:length(y_test_all)) #shuffle test set
y_test_all <- y_test_all[perm]
muhat_test_all <- muhat_test_all[perm]


alpha_set <- seq(0.05,0.2,0.025)
fdp <- matrix(0,nrow=num_test_set,ncol=length(alpha_set))
tpp <- matrix(0,nrow=num_test_set,ncol=length(alpha_set))

risk_upper <- matrix(0,nrow=num_test_set,ncol=length(alpha_set))
risk_lower <- matrix(0,nrow=num_test_set,ncol=length(alpha_set))
coverage <- matrix(0,nrow=num_test_set,ncol=length(alpha_set))
width <- matrix(0,nrow=num_test_set,ncol=length(alpha_set))

n_1 <- floor(n/2)
test_pairs <- t(combn(m,2))

pb <- txtProgressBar(min = 0, max = num_test_set, style = 3)

for(k in 1:num_test_set)
{
  muhat <- abs(muhat_cal_all[(n*(k-1)+1):(n*k)])
  Y <- y_cal_all[(n*(k-1)+1):(n*k)]
  
  n_1 <- floor(n/2)
  cal_abs_diff <- abs(muhat[1:n_1]-muhat[(n_1+1):(2*n_1)])
  cal_Y_null <- ( (muhat[1:n_1]-muhat[(n_1+1):(2*n_1)])*(Y[1:n_1]-Y[(n_1+1):(2*n_1)]) <= 0 )
  
  muhat_test <- abs(muhat_test_all[(m*(k-1)+1):(m*k)])
  Y_test <- y_test_all[(m*(k-1)+1):(m*k)]
  test_abs_diff <- apply(test_pairs,1, function(a) abs(muhat_test[a[1]] - muhat_test[a[2]]))
  
  test_null <- apply(test_pairs,1, function(a) (muhat_test[a[1]] - muhat_test[a[2]])*(Y_test[a[1]]-Y_test[a[2]]) <= 0)
  test_non_null <- 1-test_null
  
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
  
  true_rank <- rank(Y_test[order(muhat_test)])
  
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
    
    #construct uncertainty graph with e-BH procedure
    
    e_val <- ((test_abs_diff > T_rej)/(sum((cal_abs_diff>T_rej)*cal_Y_null)+1))*(n_1+1)
    
    rej_pairs <- (e_val >= eBH(e_val,alpha))
    
    R <- sum(rej_pairs)
    V <- sum(rej_pairs*test_null)
    pos <- sum(test_non_null)
    tr_pos <- sum(rej_pairs*test_non_null)
    
    fdp[k,j] <- V/max(1,R)
    tpp[k,j] <- tr_pos/max(1,pos)
    
    #prediction interval for the ranks
    
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

fdr <- colMeans(fdp)
fdr_se <- apply(fdp,2,sd)/sqrt(num_test_set)
power <- colMeans(tpp)
power_se <- apply(tpp,2,sd)/sqrt(num_test_set)

risk_upper_mean <- colMeans(risk_upper)
risk_upper_se <- apply(risk_upper,2,sd)/sqrt(num_test_set)
risk_lower_mean <- colMeans(risk_lower)
risk_lower_se <- apply(risk_lower,2,sd)/sqrt(num_test_set)

coverage_mean <- colMeans(coverage)
coverage_se <- apply(coverage,2,sd)/sqrt(num_test_set)
width_mean <- colMeans(width)
width_se <- apply(width,2,sd)/sqrt(num_test_set)


#plot
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

plot(1-alpha_set,coverage_mean,pch=16,ylim=c(0.5,1),xlab=expression(1-alpha),ylab="Mean coverage rate",cex=0.8)
arrows(1-alpha_set, coverage_mean-coverage_se, 1-alpha_set, coverage_mean+coverage_se, length=0, angle=90, code=3,lwd=2)
segments(0.8,0.8,0.95,0.95,lty=2,lwd=2)
plot(alpha_set,width_mean,pch=16,ylim=c(0,20),xlab=expression(alpha),ylab="Mean width",cex=0.8)
arrows(alpha_set, width_mean-width_se, alpha_set, width_mean+width_se, length=0, angle=90, code=3,lwd=2)


