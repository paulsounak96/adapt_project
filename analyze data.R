library(ggplot2)
library(dplyr)
library(shiny)
library(ggplot2)
library(dplyr)
library(shiny)
library(reshape2)

results <- read.csv('results.csv')
results2 <- read.csv('resultsBIG.csv')


results[results$rho ==0,] <- results
        
df <- melt(data = results, id.vars = c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
df2 <- melt(data = results2, id.vars = c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
df <- rbind(df,df2)

list <- strsplit(as.character(df$variable),"\\.")
t<- t(unlist(sapply(list, function(x) c(x[1],x[2]))))
df <- data.frame(df,t)

colnames(df)[12:13] <- c("metric","method")
agg <- aggregate(df$value,by = list(t.fdr=df$t.fdr,n=df$n,p=df$p,spr=df$spr,s=df$s,k=df$k,rho=df$rho,amplitude=df$amplitude,metric=df$metric,method=df$method),FUN = mean)
n=6000
spr=0.15
rho=0
s=0.5
metric = "power"

dat <- agg[agg["n"] ==n & agg["spr"] ==spr & agg["s"] == s & agg["rho"] == rho & agg["metric"]== metric,]

dat %>%  ggplot( aes(x=dat$t.fdr, y=dat$x, group=dat$method, color=dat$method)) +
  geom_line() + xlab("Target FDR") + ylab("Power") +ggtitle("Power") + 
  theme(plot.title = element_text(size=22,hjust = 0.5),legend.title = element_blank()) 

numtrials <- aggregate(df$value,by = list(t.fdr=df$t.fdr,n=df$n,p=df$p,spr=df$spr,s=df$s,k=df$k,rho=df$rho,amplitude=df$amplitude,metric=df$metric,method=df$method),FUN = length)
numtrials <- numtrials[numtrials$t.fdr ==0.01 & numtrials$metric=='fdr'&numtrials$method=="a",]
