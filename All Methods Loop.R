set.seed(666)
library(ggplot2)
library(doMC)
source("kada_code.R")
source("loop_functions.R")

# Problem parameters
n = 1000          # number of observations
p = 250    # number of variables
s = 0.5
amplitude = 5   # signal amplitude (for noise level = 1)
rho = 0
ko_stat = stat.glmnet_lambdasmax
#######################################################


res <- simulate_loop(n,p,s,amplitude,rho,ko_stat=ko_stat,ntrials=2,out="test.csv")

nlist <- c(1000)
ndivplist <- c(0.1,0.25)
slist <- c(0.1,0.25,0.5)
rholist <- c(0,0.1,0.2,0.3,0.4,0.5)
amplitudelist <- c(5)
ko_stat = stat.glmnet_lambdasmax
iter_list <- expand.grid(nlist,ndivplist,slist,rholist,amplitudelist)
outfile="results.csv"
ntrials = 50

for(i in 1:length(iter_list)){
  n = iter_list[i,1]    
  p = n*iter_list[i,2]     
  s = iter_list[i,3]
  rho = iter_list[i,4]
  amplitude = iter_list[i,5]  
  
  t<- try(simulate_loop(n,p,s,amplitude,rho,ko_stat=ko_stat,ntrials=ntrials,out=outfile))
}
