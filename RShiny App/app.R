library(ggplot2)
library(dplyr)
library(shiny)
library(reshape2)


results <- read.csv('results.csv')
results2 <- read.csv('resultsBIG.csv')
results3 <- read.csv('resultsMAC2.csv')
results4 <- read.csv('resultsMACnew.csv')
results5 <- read.csv('simulations.csv')

results <- results[results$rho ==0,] 
results2 <- results2[results2$rho ==0,] 
results3 <- results3[results3$rho ==0,]

df <- melt(data = results, id.vars = c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
df2 <- melt(data = results2, id.vars = c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
df3 <- melt(data = results3, id.vars = c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
df4 <- melt(data = results4, id.vars = c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
df5 <- melt(data = results5, id.vars = c('n','p','spr','s','k','rho','amplitude','ntrials','t.fdr'))
df <- rbind(df,df2,df3,df4,df5)

list <- strsplit(as.character(df$variable),"\\.")
t<- t(unlist(sapply(list, function(x) c(x[1],x[2]))))
df <- data.frame(df,t)

colnames(df)[12:13] <- c("metric","method")
agg <- aggregate(df$value,by = list(t.fdr=df$t.fdr,n=df$n,p=df$p,spr=df$spr,s=df$s,k=df$k,rho=df$rho,amplitude=df$amplitude,metric=df$metric,method=df$method),FUN = function(x) c(mean(x),length(x)))
agg$mean <- agg$x[,1]
agg$ntrials <- agg$x[,2]
agg <- agg[agg$ntrials>10,]
write.csv(agg,'aggregated_results.csv')

ui <- fluidPage(
  titlePanel("AdaPT Variations for Dependent Data"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Choose dataset parameters"),
      
      selectInput("n", 
                  label = "Number of Data Points (n)",
                  choices = unique(agg$n)),
      
      selectInput("spr", 
                  label = "Number of Variables (p/n)",
                  choices = unique(agg$spr)),
      selectInput("s", 
                  label = "Percentage Non-Noise Variables (s)",
                  choices = unique(agg$s)),
      selectInput("rho", 
                  label = "Level of correlation (rho)",
                  choices = unique(agg$rho))
      
    ),
    
    mainPanel(tabPanel("Plot",plotOutput("chartPOWER"),plotOutput('chartFDR')))
  )
)

server <- function(input, output) {
  output$chartPOWER <- renderPlot({
    par(mfrow=c(1,2))
    
    dat <- agg[agg["n"] ==input$n & agg["spr"] ==input$spr & agg["s"] == input$s & agg["rho"] == input$rho & agg["metric"]== "power",]
    
    dat %>%  ggplot( aes(x=dat$t.fdr, y=dat$mean, group=dat$method, color=dat$method)) +
      geom_line() + xlab("Target FDR") + ylab("Power") +ggtitle("Power") + 
      theme(plot.title = element_text(size=22,hjust = 0.5),legend.title = element_blank()) 
    
    
  })
  
  output$chartFDR <- renderPlot({
    par(mfrow=c(1,2))
    
    dat <- agg[agg["n"] ==input$n & agg["spr"] ==input$spr & agg["s"] == input$s & agg["rho"] == input$rho & agg["metric"]== "fdr",]
    
    dat %>%  ggplot( aes(x=dat$t.fdr, y=dat$mean, group=dat$method, color=dat$method)) +
      geom_line() + xlab("Target FDR") + ylab("False Discovery Rate") +ggtitle("False Discovery Rate") + 
      theme(plot.title = element_text(size=22,hjust = 0.5),legend.title = element_blank()) 
    
    
  })
}

shinyApp(ui = ui, server = server)

