set.seed(666)
library(ggplot2)
library(doMC)
source("kada_code.R")
source("loop_functions.R")


nlist <- c(5000,6000)
ndivplist <- c(0.1,0.15,0.25)
slist <- c(0.1,0.25,0.5)
rholist <- c(0,0.025,0.05,0.075,0.1,0.2,0.3,0.4,0.5)
amplitudelist <- c(5)
ko_stat = stat.glmnet_lambdasmax
iter_list <- expand.grid(nlist,ndivplist,slist,rholist,amplitudelist)
iter_list1 <- iter_list[order(iter_list[,4],iter_list[,1],-iter_list[,2],-iter_list[,3]),]
iter_list <- iter_list[order(iter_list[,1],iter_list[,4],-iter_list[,2],-iter_list[,3]),]

outfile="results5.csv"
ntrials = 5

for(i in 10:nrow(iter_list)){
  n = iter_list[i,1]    
  p = n*iter_list[i,2]     
  s = iter_list[i,3]
  rho = iter_list[i,4]
  amplitude = iter_list[i,5]  
  
  t<- try(simulate_loop(n,p,s,amplitude,rho,ko_stat=ko_stat,ntrials=ntrials,out=outfile))
}
simulate_loop(n,p,s,amplitude,rho,ko_stat=ko_stat,ntrials=ntrials,out=outfile)


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
