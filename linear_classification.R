linclass <- function(n,g,a,beta,sd,PI,mu,Sigma){
  ############################参数说明###############################
  ## inputs: n: number of observations(number)
  ##         g: number of groups(number)
  ##         a: intercepts of the g regression models(vector)
  ##         beta: coefficients of the g regressions(list)
  ##         sd: standard deviation of the residuals(vector)
  ##         PI: probability of each mixture component(vector)
  ##         mu: mean vectors of g X's(matrix_p*g)
  ##         Sigma: covariance matrices of g components(array_p*p*g)
  ###################################################################
  
  

  ############linear model ###############
  source('r_Gaussian_mixed.R')#generate random vectors from Gaussian mixed distribution
  source('gmm_clust1d.R') #gmm cluster by EM
  library(MASS) #use the function mvrnorn to generate multivariate normal vector
  library(dr) #using sir to reduce dimension
  
  
  p <- dim(mu)[1] 
  D <- data.frame(r_Gaussian_mixed(n=600,mu=mu,sigma = Sigma,
                                   weight = PI ,tag = TRUE)) #design matrix
  colnames(D) <- c(paste('x',1:p,sep = ""),"group")
  x <- apply(data.frame(c(1:g)), 1, 
             function(x){as.matrix(subset(D,group==x)[,-(p+1)])}) #g groups of design matrix
  
  
  ######generate Y's#########
  y <- list()
  frame <- data.frame()
  ## y_j=a_j+X*beta_j+e_j; j=1,2,...,g
  for (i in 1:g) {
    y[[i]] <- a[i]+x[[i]]%*%beta[[i]]+rnorm(dim(x[[i]])[1],mean = 0,sd=sd[i])
    frame <- rbind(frame,data.frame(x[[i]],y=y[[i]],group=i))
  }
  
  
  ################classification#################
  clust_result <- gmm_clust1d(data = frame$y, k = g,mu = sapply(y,mean),
                              sigma2 = sapply(y,var), PI = rep(1,g)/g, N = 20)
  #the best result by using the group sample mean and sample variance as priors
  
  tab <- table(frame$group,clust_result$clust.pred)
  TPM <- 1-sum(diag(tab))/sum(tab)
  
  G <- apply(data.frame(c(1:g)), 1, 
             function(x){frame[which(clust_result$clust.pred==x),]})
  #G is a list of g groups based on clust.pred
  
  
  ###############dimension reduction#############
  xnam <- paste0("x", 1:p)
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"))) #create formula like 'y~x1+x2+x3'
  sir.out <- apply(data.frame(c(1:g)), 1, function(x){dr(fmla,data = G[[x]])}) #the sir results of g groups
  #the correlation between true beta and the predicted ones
  r <- sapply(as.list(1:g),function(x){t(beta[[x]])%*%sir.out[[x]]$evectors[,1]/norm(beta[[x]],'2')})

  return(list(TPM=TPM,sir.out=sir.out,r=r))
  
}
