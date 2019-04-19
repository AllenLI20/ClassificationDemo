
#######################参数说明############################
## inputs: data, k(# of clusters), prior.mu,             ##
##         prior.Sigma, prior.PI, iter                   ##
## outputs: post.mu, post.Sigma, post.PI, maxind, maxv   ##
###########################################################

gmm_clust1d <- function(data,k,mu,sigma2,PI,N){
  
  #######################初始化######################
  
  n=length(data) #length of the vector
  
  ##################################################
  
  if (N>0) {
    
    for (times in 1:N) {
      
      para <- cbind(PI,mu,sigma2) #put all the parameters in one dataframe(matrix)
      
      #dataframe to put the posteriors
      gamma <- apply(para, 1, function(x){x[1]*dnorm(data,x[2],sqrt(x[3]))})
      gamma <- gamma/apply(gamma, 1, sum)
      colnames(gamma) <- paste('posterior',1:k)
      maxind <- apply(gamma, 1, which.max)
      maxv <- apply(gamma,1,max)
      
      for (j in 1:k) {
        mu[j] <- sum(gamma[,j]*data)/sum(gamma[,j])
        sigma2[j] <- sum(gamma[,j]*(data-mu[j])^2)/sum(gamma[,j])
        PI[j] <- sum(gamma[,j])/n
      }
      
    }
    
    
    rtn <- list(mu=mu,sigma2=sigma2,PI=PI,clust.pred=maxind,posterior=maxv)
    return(rtn)
    
    
  }
  else return()
  
}
