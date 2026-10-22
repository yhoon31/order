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

##--------------------------------------------------------

#load data
y_cal_all <- read.table("y_cal.txt")$V1
y_test_all <- read.table("y_test.txt")$V1
muhat_cal_all <- read.table("muhat_cal.txt")$V1
muhat_test_all <- read.table("muhat_test.txt")$V1


##FDR control (main paper)
n <- 200
m <- 20
num_test_set <- 80

set.seed(1234)
perm_cal <- sample(1:length(y_cal_all))
y_cal_all <- y_cal_all[perm_cal]
muhat_cal_all <- muhat_cal_all[perm_cal]
perm <- sample(1:length(y_test_all))
y_test_all <- y_test_all[perm]
muhat_test_all <- muhat_test_all[perm]

alpha_set <- seq(0.05,0.2,0.025)
n_alpha <- length(alpha_set)

fdp <- matrix(0,nrow=num_test_set,ncol=n_alpha)
tpp <- matrix(0,nrow=num_test_set,ncol=n_alpha)
fdp_U <- matrix(0,nrow=num_test_set,ncol=n_alpha)
tpp_U <- matrix(0,nrow=num_test_set,ncol=n_alpha)
fdp_p <- matrix(0,nrow=num_test_set,ncol=n_alpha)
tpp_p <- matrix(0,nrow=num_test_set,ncol=n_alpha)

risk_upper <- matrix(0,nrow=num_test_set,ncol=n_alpha)
risk_lower <- matrix(0,nrow=num_test_set,ncol=n_alpha)
coverage <- matrix(0,nrow=num_test_set,ncol=n_alpha)
width <- matrix(0,nrow=num_test_set,ncol=n_alpha)

risk_upper_U <- matrix(0,nrow=num_test_set,ncol=n_alpha)
risk_lower_U <- matrix(0,nrow=num_test_set,ncol=n_alpha)
coverage_U <- matrix(0,nrow=num_test_set,ncol=n_alpha)
width_U <- matrix(0,nrow=num_test_set,ncol=n_alpha)

risk_upper_p <- matrix(0,nrow=num_test_set,ncol=n_alpha)
risk_lower_p <- matrix(0,nrow=num_test_set,ncol=n_alpha)
coverage_p <- matrix(0,nrow=num_test_set,ncol=n_alpha)
width_p <- matrix(0,nrow=num_test_set,ncol=n_alpha)

n_1 <- floor(n/2)
test_pairs <- t(combn(m,2))
n_pairs <- nrow(test_pairs)
pair_1 <- test_pairs[,1]
pair_2 <- test_pairs[,2]
pair_factor <- choose(m-2,2)/(n_1+1)
pairs_by_index <- lapply(1:m,function(a) which(pair_1==a | pair_2==a))

pair_index <- matrix(NA_integer_,m,m)
pair_index[cbind(pair_1,pair_2)] <- seq_along(pair_1)
pair_index[cbind(pair_2,pair_1)] <- seq_along(pair_1)

pb <- txtProgressBar(min=0,max=num_test_set,style=3)

for(k in 1:num_test_set)
{
  muhat <- abs(muhat_cal_all[(n*(k-1)+1):(n*k)])
  Y <- y_cal_all[(n*(k-1)+1):(n*k)]
  
  cal_diff <- muhat[1:n_1]-muhat[(n_1+1):(2*n_1)]
  cal_abs_diff <- abs(cal_diff)
  cal_Y_null <- cal_diff*(Y[1:n_1]-Y[(n_1+1):(2*n_1)])<=0
  cal_null_sorted <- sort(cal_abs_diff[cal_Y_null])
  n_cal_null <- length(cal_null_sorted)
  
  muhat_test <- abs(muhat_test_all[(m*(k-1)+1):(m*k)])
  muhat_rank <- rank(muhat_test + runif(m, -1e-10, 1e-10))
  
  idx_list <- vector("list",n_pairs)
  for(r in seq_len(n_pairs))
  {
    lower <- which(muhat_rank <= min(muhat_rank[pair_1[r]],muhat_rank[pair_2[r]]))
    upper <- which(muhat_rank >= max(muhat_rank[pair_1[r]],muhat_rank[pair_2[r]]))
    idx_list[[r]] <- as.vector(pair_index[lower,upper,drop=FALSE])
  }
  
  Y_test <- y_test_all[(m*(k-1)+1):(m*k)]
  
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
    risk_upper_j <- sum(pmax(true_rank-rank_interval[,2],0))/max(sum(m-rank_interval[,2]),1)
    risk_lower_j <- sum(pmax(rank_interval[,1]-true_rank,0))/max(sum(rank_interval[,1]-1),1)
    coverage_j <- mean(rank_interval[,1]<=true_rank & true_rank<=rank_interval[,2])
    width_j <- mean(rank_interval[,2]-rank_interval[,1]+1)
    c(risk_upper_j,risk_lower_j,coverage_j,width_j)
  }
  
  for(j in seq_len(n_alpha))
  {
    alpha <- alpha_set[j]
    T_rej <- T_rej_mat[,j]
    
    finite_T <- is.finite(T_rej)
    e_denom <- rep(1,n_pairs)
    e_denom[finite_T] <- n_cal_null-findInterval(T_rej[finite_T],cal_null_sorted)+1
    e_val <- as.numeric(test_abs_diff>T_rej)*(n_1+1)/e_denom
    
    e_val_tilde <- sapply(idx_list,function(idx) min(e_val[idx]))
    
    U <- runif(1)
    e_val_tilde_U <- e_val_tilde/U
    
    rej_pairs <- e_val_tilde>=eBH(e_val_tilde,alpha)
    rej_pairs_U <- e_val_tilde_U>=eBH(e_val_tilde_U,alpha)
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
    
    fdp[k,j] <- V/max(1,R)
    tpp[k,j] <- tr_pos/max(1,pos)
    fdp_U[k,j] <- V_U/max(1,R_U)
    tpp_U[k,j] <- tr_pos_U/max(1,pos)
    fdp_p[k,j] <- V_p/max(1,R_p)
    tpp_p[k,j] <- tr_pos_p/max(1,pos)
    
    rank_res <- get_rank_result(rej_pairs)
    rank_res_U <- get_rank_result(rej_pairs_U)
    rank_res_p <- get_rank_result(rej_pairs_p)
    
    risk_upper[k,j] <- rank_res[1]
    risk_lower[k,j] <- rank_res[2]
    coverage[k,j] <- rank_res[3]
    width[k,j] <- rank_res[4]
    
    risk_upper_U[k,j] <- rank_res_U[1]
    risk_lower_U[k,j] <- rank_res_U[2]
    coverage_U[k,j] <- rank_res_U[3]
    width_U[k,j] <- rank_res_U[4]
    
    risk_upper_p[k,j] <- rank_res_p[1]
    risk_lower_p[k,j] <- rank_res_p[2]
    coverage_p[k,j] <- rank_res_p[3]
    width_p[k,j] <- rank_res_p[4]
  }
  
  setTxtProgressBar(pb,k)
}

close(pb)

fdr <- colMeans(fdp)
fdr_se <- apply(fdp,2,sd)/sqrt(num_test_set)
power <- colMeans(tpp)
power_se <- apply(tpp,2,sd)/sqrt(num_test_set)

fdr_U <- colMeans(fdp_U)
fdr_U_se <- apply(fdp_U,2,sd)/sqrt(num_test_set)
power_U <- colMeans(tpp_U)
power_U_se <- apply(tpp_U,2,sd)/sqrt(num_test_set)

fdr_p <- colMeans(fdp_p)
fdr_p_se <- apply(fdp_p,2,sd)/sqrt(num_test_set)
power_p <- colMeans(tpp_p)
power_p_se <- apply(tpp_p,2,sd)/sqrt(num_test_set)

risk_upper_mean <- colMeans(risk_upper)
risk_upper_se <- apply(risk_upper,2,sd)/sqrt(num_test_set)
risk_lower_mean <- colMeans(risk_lower)
risk_lower_se <- apply(risk_lower,2,sd)/sqrt(num_test_set)
coverage_mean <- colMeans(coverage)
coverage_se <- apply(coverage,2,sd)/sqrt(num_test_set)
width_mean <- colMeans(width)
width_se <- apply(width,2,sd)/sqrt(num_test_set)

risk_upper_U_mean <- colMeans(risk_upper_U)
risk_upper_U_se <- apply(risk_upper_U,2,sd)/sqrt(num_test_set)
risk_lower_U_mean <- colMeans(risk_lower_U)
risk_lower_U_se <- apply(risk_lower_U,2,sd)/sqrt(num_test_set)
coverage_U_mean <- colMeans(coverage_U)
coverage_U_se <- apply(coverage_U,2,sd)/sqrt(num_test_set)
width_U_mean <- colMeans(width_U)
width_U_se <- apply(width_U,2,sd)/sqrt(num_test_set)

risk_upper_p_mean <- colMeans(risk_upper_p)
risk_upper_p_se <- apply(risk_upper_p,2,sd)/sqrt(num_test_set)
risk_lower_p_mean <- colMeans(risk_lower_p)
risk_lower_p_se <- apply(risk_lower_p,2,sd)/sqrt(num_test_set)
coverage_p_mean <- colMeans(coverage_p)
coverage_p_se <- apply(coverage_p,2,sd)/sqrt(num_test_set)
width_p_mean <- colMeans(width_p)
width_p_se <- apply(width_p,2,sd)/sqrt(num_test_set)

res <- cbind(fdr,fdr_se,fdr_U,fdr_U_se,fdr_p,fdr_p_se,
             power,power_se,power_U,power_U_se,power_p,power_p_se,
             risk_upper_mean,risk_upper_se,risk_upper_U_mean,risk_upper_U_se,risk_upper_p_mean,risk_upper_p_se,
             risk_lower_mean,risk_lower_se,risk_lower_U_mean,risk_lower_U_se,risk_lower_p_mean,risk_lower_p_se,
             coverage_mean,coverage_se,coverage_U_mean,coverage_U_se,coverage_p_mean,coverage_p_se,
             width_mean,width_se,width_U_mean,width_U_se,width_p_mean,width_p_se)


col1 <- "#0072B2"
col3 <- "#B2182B"
col2 <- "#56B4E9"
pch_set <- c(16,16,4)

par(mfrow=c(3,2),cex.lab=2,cex.axis=1.7,mar=c(5,6,2,0),xpd=NA,oma=c(0,0,0,4),mgp=c(3,1,0))

matplot(alpha_set,res[,c(1,3,5)],pch=pch_set,ylim=c(0,0.22),xlab=expression(alpha),ylab="FDR",cex=0.8,col=c(col1,col2,col3))
arrows(alpha_set,res[,1]-res[,2],alpha_set,res[,1]+res[,2],length=0,angle=90,code=3,lwd=2,col=col1)
arrows(alpha_set,res[,3]-res[,4],alpha_set,res[,3]+res[,4],length=0,angle=90,code=3,lwd=2,col=col2)
arrows(alpha_set,res[,5]-res[,6],alpha_set,res[,5]+res[,6],length=0,angle=90,code=3,lwd=2,col=col3)
segments(0.05,0.05,0.2,0.2,lty=2,lwd=2)
legend("topleft",legend=c("e-BH","U-eBH","p-BH"),pch=pch_set,col=c(col1,col2,col3),bty="n")

matplot(alpha_set,res[,c(7,9,11)],pch=pch_set,ylim=c(0,1),xlab=expression(alpha),ylab="Power",cex=0.8,col=c(col1,col2,col3))
arrows(alpha_set,res[,7]-res[,8],alpha_set,res[,7]+res[,8],length=0,angle=90,code=3,lwd=2,col=col1)
arrows(alpha_set,res[,9]-res[,10],alpha_set,res[,9]+res[,10],length=0,angle=90,code=3,lwd=2,col=col2)
arrows(alpha_set,res[,11]-res[,12],alpha_set,res[,11]+res[,12],length=0,angle=90,code=3,lwd=2,col=col3)

matplot(alpha_set,res[,c(13,15,17)],pch=pch_set,ylim=c(0,0.22),xlab=expression(alpha),ylab="Risk+",cex=0.8,col=c(col1,col2,col3))
arrows(alpha_set,res[,13]-res[,14],alpha_set,res[,13]+res[,14],length=0,angle=90,code=3,lwd=2,col=col1)
arrows(alpha_set,res[,15]-res[,16],alpha_set,res[,15]+res[,16],length=0,angle=90,code=3,lwd=2,col=col2)
arrows(alpha_set,res[,17]-res[,18],alpha_set,res[,17]+res[,18],length=0,angle=90,code=3,lwd=2,col=col3)
segments(0.05,0.05,0.2,0.2,lty=2,lwd=2)

matplot(alpha_set,res[,c(19,21,23)],pch=pch_set,ylim=c(0,0.22),xlab=expression(alpha),ylab="Risk-",cex=0.8,col=c(col1,col2,col3))
arrows(alpha_set,res[,19]-res[,20],alpha_set,res[,19]+res[,20],length=0,angle=90,code=3,lwd=2,col=col1)
arrows(alpha_set,res[,21]-res[,22],alpha_set,res[,21]+res[,22],length=0,angle=90,code=3,lwd=2,col=col2)
arrows(alpha_set,res[,23]-res[,24],alpha_set,res[,23]+res[,24],length=0,angle=90,code=3,lwd=2,col=col3)
segments(0.05,0.05,0.2,0.2,lty=2,lwd=2)

matplot(1-alpha_set,res[,c(25,27,29)],pch=pch_set,ylim=c(0.5,1),xlab=expression(1-alpha),ylab="Mean coverage rate",cex=0.8,col=c(col1,col2,col3))
arrows(1-alpha_set,res[,25]-res[,26],1-alpha_set,res[,25]+res[,26],length=0,angle=90,code=3,lwd=2,col=col1)
arrows(1-alpha_set,res[,27]-res[,28],1-alpha_set,res[,27]+res[,28],length=0,angle=90,code=3,lwd=2,col=col2)
arrows(1-alpha_set,res[,29]-res[,30],1-alpha_set,res[,29]+res[,30],length=0,angle=90,code=3,lwd=2,col=col3)
segments(0.8,0.8,0.95,0.95,lty=2,lwd=2)

matplot(alpha_set,res[,c(31,33,35)],pch=pch_set,ylim=c(0,m),xlab=expression(alpha),ylab="Mean width",cex=0.8,col=c(col1,col2,col3))
arrows(alpha_set,res[,31]-res[,32],alpha_set,res[,31]+res[,32],length=0,angle=90,code=3,lwd=2,col=col1)
arrows(alpha_set,res[,33]-res[,34],alpha_set,res[,33]+res[,34],length=0,angle=90,code=3,lwd=2,col=col2)
arrows(alpha_set,res[,35]-res[,36],alpha_set,res[,35]+res[,36],length=0,angle=90,code=3,lwd=2,col=col3)



##FWER control (Appendix)
y_cal <- read.table("y_cal.txt")$V1
y_test_all <- read.table("y_test.txt")$V1
muhat <- read.table("muhat_cal.txt")$V1
muhat_test_all <- read.table("muhat_test.txt")$V1

set.seed(1)
perm <- sample(1:length(y_test_all)) #shuffle test set
y_test_all <- y_test_all[perm]
muhat_test_all <- muhat_test_all[perm]

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
