set.seed(666)
library(ggplot2)
library(doMC)
source("kada_code.R")
source("loop_functions.R")


nlist <- c(1000,5000,6000)
ndivplist <- c(0.1,0.15,0.25)
slist <- c(0.1,0.25,0.5)
rholist <- c(0,0.025,0.05,0.075,0.1,0.2,0.3,0.4,0.5)
amplitudelist <- c(5)
ko_stat = stat.glmnet_lambdasmax
iter_list <- expand.grid(nlist,ndivplist,slist,rholist,amplitudelist)
iter_list <- iter_list[order(iter_list[,4],iter_list[,1],-iter_list[,2],-iter_list[,3]),]

outfile="results2.csv"
ntrials = 50

for(i in 6:length(iter_list)){
  n = iter_list[i,1]    
  p = n*iter_list[i,2]     
  s = iter_list[i,3]
  rho = iter_list[i,4]
  amplitude = iter_list[i,5]  
  
  t<- try(simulate_loop(n,p,s,amplitude,rho,ko_stat=ko_stat,ntrials=ntrials,out=outfile))
}
