########################## FRANKENSTEIN (without double dipping) ###########################

# Problem parameters
n = 5000          # number of observations
p = 1000           # number of variables
k = 100           # number of variables with nonzero coefficients
amplitude = 4.5   # signal amplitude (for noise level = 1)

# Generate the variables from a multivariate normal distribution
mu = rep(0,p)
rho = 0.25
Sigma = toeplitz(rho^(0:(p-1)))
X = matrix(rnorm(n*p),n) %*% chol(Sigma)

# Generate the response from a linear model
nonzero = sample(p, k)
beta = amplitude * (1:p %in% nonzero) / sqrt(n)
y.sample = function(X) X %*% beta + rnorm(n)
y = y.sample(X)



# Helper functions
ROC_point <- function(result,alpha){
  if (class(result) == "knockoff.result"){
    t <- knockoff.threshold(result$statistic,fdr=alpha)
    c <- c(alpha,power(which((result$statistic>t) == TRUE)),fdr(which((result$statistic>t) == TRUE)))
  }
  if (class(result) == "adapt"){
    c = c(alpha, power(as.numeric(result$rejs[[100*alpha]])), fdr(as.numeric(result$rejs[[100*alpha]])))
  }
  
  return(c)
}
ROC_curve <- function(result,alphas=seq(0.01,1,by=0.01)){
  t(sapply(alphas,function(x) ROC_point(result,x)))
}
fdr = function(selected) sum(beta[selected] == 0) / max(1, length(selected))
power = function(selected) sum(beta[selected] != 0) / k



# Partition the dataset
samp = sample(n, n/2)
X1 = X[samp, ]
y1 = y[samp]
X2 = X[-samp, ]
y2 = y[-samp]



# Knockoffs first
result_knockoff = knockoff.filter(X, y, knockoffs=create.fixed, statistic=stat.glmnet_lambdasmax)

result_temp = knockoff.filter(X2, y2, knockoffs=create.fixed,
                              statistic=stat.glmnet_lambdasmax)



# Then AdaPT
formulas <- paste0("ns(x, df = ", 6:10, ")")
models <- lapply(formulas, function(formula){
  piargs <- muargs <- list(formula = formula)
  gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
})

pvals_adapt = summary(lm(y~X))[["coefficients"]][,"Pr(>|t|)"][-1]
x_adapt = data.frame(x = as.numeric(summary(lm(y~X))[["coefficients"]][,"Std. Error"][-1]))
result_adapt <- adapt(x = x_adapt, pvals = pvals_adapt, models = models, dist = beta_family(), nfits = 10)

pvals_temp = summary(lm(y1~X1))[["coefficients"]][,"Pr(>|t|)"][-1]
x_temp <- data.frame(x = as.numeric(result_temp$statistic))
result_franken <- adapt(x = x_temp, pvals = pvals_temp, models = models, dist = beta_family(), nfits = 10)



# ROC Curves and plotting
roc_knockoff = ROC_curve(result_knockoff)
roc_adapt = ROC_curve(result_adapt)
roc_franken = ROC_curve(result_franken)

par(mfrow=c(1,2))
plot(roc_knockoff[,1],roc_knockoff[,2],type="l",main="Power, rho = 0.25",xlab="Target FDR Level",ylab="Power",ylim=c(0,1))
legend(0.1,0.2,c("Knockoff","AdaPT with Knockoffs as side info", "AdaPT"),lty=c(1,2,3))
lines(roc_franken[,1],roc_franken[,2],type="l",lty=2)
lines(roc_adapt[,1],roc_adapt[,2],type="l",lty=3)

plot(roc_knockoff[,1],roc_knockoff[,3],type="l",main="Level, rho = 0.25",xlab="Target FDR Level",ylab="Level",ylim=c(0,1))
lines(roc_franken[,1],roc_franken[,3],type="l",lty=2)
lines(roc_adapt[,1],roc_adapt[,3],type="l",lty=3)