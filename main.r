set.seed(666)
library(ggplot2)
library(doMC)
library(knockoff)
library(kableExtra)



nlist <- c(6000)
ndivplist <- c(0.15,0.2)
slist <- c(0.1)
rholist <- c(0,0.05,0.1,0.15,0.2)
amplitudelist <- c(5)
iter_list <- expand.grid(nlist,ndivplist,slist,rholist,amplitudelist)
ntrials = 10
outfile = "simulations.csv"

for(mnum in 1:nrow(iter_list)){
  
  n = iter_list[mnum,1]    
  spr = iter_list[mnum,2]  
  p = n*spr    
  s = iter_list[mnum,3]
  rho = iter_list[mnum,4]
  amplitude = iter_list[mnum,5]  
  k = p * s   
  
  setup <- data.frame(n,p,spr,s,k,rho,amplitude,ntrials,t.fdr=1:20 * 0.01)
  
  dpower.kadapt <- setup
  dpower.k <- setup
  dpower.a <- setup
  dpower.sounak <- setup
  dpower.james <- setup
  
  dfdr.kadapt <- setup
  dfdr.k <- setup
  dfdr.a <- setup
  dfdr.sounak <- setup
  dfdr.james <- setup
  
  for(t in 1:ntrials) {
    #---------------------------------------------------------------
    # Adaptive P-Value Thresholding For Multiple Testing With Side Information
    #---------------------------------------------------------------
    # Problem parameters
    # Generate the variables from a multivariate normal distribution
    start <- Sys.time()
    source("kadapt.r")
    
    mu = rep(0,p)
    Sigma = toeplitz(rho^(0:(p-1)))
    X = matrix(rnorm(n*p),n) %*% chol(Sigma)
    
    # Generate the response from a linear model
    nonzero = sample(p, k)
    beta = amplitude * (1:p %in% nonzero) / sqrt(n)
    y.sample = function(X) X %*% beta + rnorm(n)
    y = y.sample(X)
    
    
    create_knockoffs = function(X) {
      create.second_order(X, shrink=T)
    }
    
    knockoffs = function(X) create.second_order(X, method='equi')
    
    result <- knockoff.filter(X, y, knockoffs = create_knockoffs)
    
    ti_X <- cbind(result$X, result$Xk)
    dim(ti_X)
    
    
    
    mod.1 <- lm(y ~ ti_X)
    p_tot <- summary(mod.1)$coefficients[,4]  
    
    #p_0 <- p[1:(length(p)/2)]
    #p_t <- p[(length(p)/2 + 1):length(p)]
    
    mod.2 <- lm(y ~ X)
    kable(head(summary(mod.2)$coefficients),"markdown")
    
    p_0 <- summary(mod.2)$coefficients[,4][-1]  
    
    mod.3 <- lm(y ~ result$Xk)
    p_t <- summary(mod.3)$coefficients[,4][-1]
    
    
    
    # side cov
    x <- summary(mod.2)$coefficients[,2][-1]
    x <- data.frame(x = as.numeric(x))
    
    
    # Generate models for function adapt
    library("splines")
    formulas <- paste0("ns(x, df = ", 6:10, ")")
    models <- lapply(formulas, function(formula){
      piargs <- muargs <- list(formula = formula)
      gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
    })
    
    
    res <- Kadapt(x = x, pvals = p_0, kpvals = p_t, models = models, dist = beta_family(), nfits = 10)
    
    
    
    
    rej.05 <- res$rejs[[5]]
    sum(rej.05 %in% nonzero) / length(rej.05)
    length(rej.05) / length(nonzero)
    
    t.fdr <- numeric()
    power.kadapt <- numeric()
    fdr.kadapt <- numeric()
    
    for(i in 1:20){
      rej <- res$rejs[[i]]
      power.kadapt <- c(power.kadapt, sum(rej %in% nonzero) / length(nonzero))
      fdr.kadapt <- c(fdr.kadapt,sum(beta[rej] == 0) / max(1, length(rej)))
    }
    
    
    library(adaptMT)
    
    t.fdr.k <- numeric()
    power.k <- numeric()
    fdr.k <- numeric()
    
    X_k = create.gaussian(X, mu, Sigma)
    W = stat.glmnet_coefdiff(X, X_k, y)
    
    
    for(i in 1:20){
      thres = knockoff.threshold(W, fdr=i * 0.01) 
      rej <- which(W >= thres)
      power.k <- c(power.k, sum(rej %in% nonzero) / length(nonzero))
      fdr.k <- c(fdr.k,sum(beta[rej] == 0) / max(1, length(rej)))
    }
    
    
    
    library(adaptMT)
    
    t.fdr.a <- numeric()
    power.a <- numeric()
    fdr.a <- numeric()
    
    
    resa <- adapt(x = x, pvals = p_0, models = models, dist = beta_family(), nfits = 10)
    
    
    
    for(i in 1:20){
      rej <- resa$rejs[[i]]
      power.a <- c(power.a, sum(rej %in% nonzero) / length(nonzero))
      fdr.a <- c(fdr.a,sum(beta[rej] == 0) / max(1, length(rej)))
    }
    
    #####################################################################################
    #Sounak Method
    ######################################################################################
    samp = sample(n, n/2)
    X1 = X[samp,]
    y1 = y[samp]
    X2 = X[-samp,]
    y2 = y[-samp]
    
    reg_split = summary(lm(y1~X1))
    result_split = knockoff.filter(X2, y2, knockoffs=create.fixed)
    
    
    pvals_split = reg_split[["coefficients"]][,"Pr(>|t|)"][-1]
    x = data.frame(x=result_split$statistic)
    res.sounak <- adapt_glm(x = x, pvals = pvals_split, pi_formulas = formulas,
                            mu_formulas = formulas,  dist = beta_family(), nfits = 10)
    
    
    t.fdr.sounak <- numeric()
    power.sounak <- numeric()
    fdr.sounak <- numeric()
    
    for(i in 1:20){
      rej <- res.sounak$rejs[[i]]
      power.sounak <- c(power.sounak, sum(rej %in% nonzero) / length(nonzero))
      fdr.sounak <- c(fdr.sounak,sum(beta[rej] == 0) / max(1, length(rej)))
    }
    
    #####################################################################################
    #James Method
    ######################################################################################
    result_james = knockoff.filter(X, y, knockoffs=create.fixed)
    
    d<- t(result_james$X)%*%result_james$X - t(result_james$Xk) %*%result_james$X
    z <- (t(result_james$X) %*% y- t(result_james$Xk) %*% y) /(sqrt(2*diag(d)))
    pvals <- 2*pmin(pnorm(z), 1- pnorm(z))
    x <- data.frame(x= (t(result_james$X) %*% y+ t(result_james$Xk) %*% y)) 
    stat_cust <- abs(t(result_james$X) %*% y)- abs(t(result_james$Xk) %*% y)
    
    res.james <- adapt_glm(x = x, pvals = pvals, pi_formulas = formulas,
                           mu_formulas = formulas, dist = beta_family(), nfits = 10)
    
    t.fdr.james <- numeric()
    power.james <- numeric()
    fdr.james <- numeric()
    
    for(i in 1:20){
      rej <- res.james$rejs[[i]]
      power.james <- c(power.james, sum(rej %in% nonzero) / length(nonzero))
      fdr.james <- c(fdr.james,sum(beta[rej] == 0) / max(1, length(rej)))
    }
    
    
    dfdr.kadapt <- data.frame(dfdr.kadapt, fdr.kadapt)
    dfdr.k <- data.frame(dfdr.k, fdr.k)
    dfdr.a <- data.frame( dfdr.a , fdr.a)
    dfdr.james <- data.frame(dfdr.james, fdr.james)
    dfdr.sounak <- data.frame( dfdr.sounak, fdr.sounak)
    
    dpower.kadapt <- data.frame(dpower.kadapt, power.kadapt)
    dpower.k <- data.frame(dpower.k, power.k)
    dpower.a <- data.frame( dpower.a , power.a)
    dpower.james <- data.frame(dpower.james, power.james)
    dpower.sounak <- data.frame( dpower.sounak, power.sounak)
    
    end <- Sys.time()
    print(end - start)
  }
  
  results_combined <- merge(dfdr.kadapt,dfdr.k,by=c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
  results_combined <- merge(results_combined,dfdr.a,,by=c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
  results_combined <- merge(results_combined,dfdr.james,,by=c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
  results_combined <- merge(results_combined,dfdr.sounak,,by=c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
  results_combined <- merge(results_combined,dpower.kadapt,,by=c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
  results_combined <- merge(results_combined,dpower.k,,by=c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
  results_combined <- merge(results_combined,dpower.a,,by=c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
  results_combined <- merge(results_combined,dpower.james,,by=c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
  results_combined <- merge(results_combined,dpower.sounak,,by=c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
  
  if(file.exists(outfile)) {
    write.table(results_combined, outfile, sep = ",",append = T,row.names=F,col.names=F)
  }
  if(!file.exists(outfile)) {
    write.table(results_combined, outfile, sep = ",",row.names=F,col.names=T)
  }
  
}

