


fdr = function(selected) sum(beta[selected] == 0) / max(1, length(selected))

power = function(selected) sum(beta[selected] != 0) / k


KO_ROC_point <- function(stat,alpha){
  t <- knockoff.threshold(stat,fdr=alpha)
  c <- c(alpha,power(which((stat>t) == TRUE)),fdr(which((stat>t) == TRUE)))
  return(c)
}


KO_ROC_curve <- function(stat,alphas=seq(0.01,1,by=0.01)){
  t(sapply(alphas,function(x) KO_ROC_point(stat,x)))
}


adapt_ROC_curve <- function(rejection_list,alphas){
  n <- length(rejection_list)
  for(i in 1:n){
    row <- c(alphas[i],power(rejection_list[[i]]),fdr(rejection_list[[i]]))
    if(i ==1){
      results <- row
    }
    else{
      results <- rbind(results,row)
    }
  }
  rownames(results) <- 1:n
  return(results)
}


simulate <- function(n,p,s,amplitude,rho,ko_stat=stat.glmnet_lambdasmax){
  set.seed(sample(1:100,1))
  k = p * s           # number of variables with nonzero coefficients
  
  # Generate the variables from a multivariate normal distribution
  mu = rep(0,p)
  Sigma = toeplitz(rho^(0:(p-1)))
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  
  # Generate the response from a linear model
  nonzero = sample(p, k)
  beta = amplitude * (1:p %in% nonzero) / sqrt(n)
  y.sample = function(X) X %*% beta + rnorm(n)
  y = y.sample(X)
  
  # Partition the dataset
  samp = sample(n, n/2)
  X1 = X[samp, ]
  y1 = y[samp]
  X2 = X[-samp, ]
  y2 = y[-samp]
  
  result = knockoff.filter(X, y, knockoffs=create.fixed,statistic=ko_stat)
  result_split = knockoff.filter(X2, y2, knockoffs=create.fixed,statistic=ko_stat)
  reg = summary(lm(y~X))
  reg_split = summary(lm(y1~X1))
  reg_ko = summary(lm(y~result$Xk))
  
  
  formulas <- paste0("ns(x, df = ", 6:10, ")")
  models <- lapply(formulas, function(formula){
    piargs <- muargs <- list(formula = formula)
    gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
  })
  
  
  ###########Standard AdaPT###########################
  library(adaptMT)
  alphas = c(seq(0.001,0.009,by=0.001),seq(0.01,1,by=0.01))
  pvals = reg[["coefficients"]][,"Pr(>|t|)"][-1]
  x <- data.frame(x=reg$coefficients[,2][-1])
  
  
  res1 <- adapt_glm(x = x, pvals = pvals, pi_formulas = formulas,
                    mu_formulas = formulas,  dist = beta_family(), nfits = 10,alphas=alphas)
  
  ###########Method 6.2 (James)############################
  d<- t(result$X)%*%result$X - t(result$Xk) %*%result$X
  z <- (t(result$X) %*% y- t(result$Xk) %*% y) /(sqrt(2*diag(d)))
  pvals <- 2*pmin(pnorm(z), 1- pnorm(z))
  x <- data.frame(x= (t(result$X) %*% y+ t(result$Xk) %*% y)) 
  stat_cust <- abs(t(result$X) %*% y)- abs(t(result$Xk) %*% y)
  
  res2 <- adapt_glm(x = x, pvals = pvals, pi_formulas = formulas,
                    mu_formulas = formulas, dist = beta_family(), nfits = 10,alphas=alphas)
  
  ###########Sounak Method############################
  pvals_split = reg_split[["coefficients"]][,"Pr(>|t|)"][-1]
  x = data.frame(x=result_split$statistic)
  res3 <- adapt_glm(x = x, pvals = pvals_split, pi_formulas = formulas,
                    mu_formulas = formulas,  dist = beta_family(), nfits = 10,alphas=alphas)
  
  ###########Yuhan Method############################
  source("kada_code.R")
  x <- data.frame(x=reg$coefficients[,2][-1])
  pvals_ko = reg_ko[["coefficients"]][,"Pr(>|t|)"][-1]
  res4 <- Kadapt(x = x, pvals = pvals, kpvals = pvals_ko, models = models, dist = beta_family(), nfits = 10)
  
  results<- cbind(KO_ROC_curve(result$statistic,alphas=alphas),"Knockoff")
  results<- rbind(results,cbind(KO_ROC_curve(stat_cust,alphas=alphas),"Knockoff (alt 1)"))
  results<-rbind(results,cbind(adapt_ROC_curve(res1$rejs,alphas=alphas),"AdaPT"))
  results<-rbind(results,cbind(adapt_ROC_curve(res2$rejs,alphas=alphas),"AdaPT (alt 1)"))
  results<-rbind(results,cbind(adapt_ROC_curve(res3$rejs,alphas=alphas),"AdaPT (alt 2)"))
  results<-rbind(results,cbind(adapt_ROC_curve(res4$rejs,alphas=alphas),"KadaPT"))
  results <-data.frame(results)
  colnames(results) <- c("Target FDR","Power","FDR","Method")
  return(results)
}


simulate_loop <- function(n,p,s,amplitude,rho,ko_stat=stat.glmnet_lambdasmax,ntrials=100,out="out.csv"){
  for(i in 1:ntrials){
    dat <- simulate(n,p,s,amplitude,rho,ko_stat=stat.glmnet_lambdasmax)
    if(class(dat) == "try-error"){
      dat <- c(rep(NA,6))
    }
    
    res <- cbind(dat,n,p,s,amplitude,rho,i,ntrials)
    
    if(file.exists(out)){
      write.table(res, out, sep = ",",append = T,row.names=F,col.names=F)
    }
    else{
      write.table(res, out, sep = ",",row.names=F,col.names=T)
    }
    
    if(i==1){
      results_compiled <- res
    }
    else {
      results_compiled <- rbind(results_compiled,res)
    }
  }
  return(results_compiled)
}

