library(ggplot2)
library(dplyr)
library(shiny)

results <- read.csv('results.csv')
agg <- aggregate(cbind(results$Power,results$FDR),by = list(results$Target.FDR,results$Method,results$n,results$p,results$s,results$amplitude,results$rho),FUN=mean)

colnames(agg) <- c("Target FDR","Method","n","p","s","Amplitude","rho","Power","FDR")
agg$pdivn <- agg['p']/agg['n']


ui <- fluidPage(
  titlePanel("AdaPT Variations for Dependent Data"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Choose dataset parameters"),
      
      selectInput("var1", 
                  label = "Number of Data Points (n)",
                  choices = unique(agg$n)),
      
      selectInput("var2", 
                  label = "Number of Variables (p/n)",
                  choices = unique(agg$pdivn)),
      selectInput("var3", 
                  label = "Percentage Non-Noise Variables (s)",
                  choices = unique(agg$s)),
      selectInput("var4", 
                  label = "Level of correlation (rho)",
                  choices = unique(agg$rho))
      
    ),
    
    mainPanel(tabPanel("Plot",plotOutput("chartPOWER"),plotOutput('chartFDR')))
  )
)

server <- function(input, output) {
  output$chartPOWER <- renderPlot({
    par(mfrow=c(1,2))
    
    dat <- agg[agg["n"] ==input$var1 & agg["pdivn"] ==input$var2 & agg["s"] == input$var3 & agg["rho"] == input$var4,]
    
    dat %>%   ggplot( aes(x=dat$`Target FDR`, y=dat$`Power`, group=dat$`Method`, color=dat$`Method`)) +
      geom_line() + xlab("Target FDR") + ylab("Power") +ggtitle("Power") + 
      theme(plot.title = element_text(size=22,hjust = 0.5),legend.title = element_blank()) 
    
    
  
    
  })
  
  output$chartFDR <- renderPlot({
    par(mfrow=c(1,2))
    
    dat <- agg[agg["n"] ==input$var1 & agg["pdivn"] ==input$var2 & agg["s"] == input$var3 & agg["rho"] == input$var4,]
    
    dat %>%   ggplot( aes(x=dat$`Target FDR`, y=dat$`FDR`, group=dat$`Method`, color=dat$`Method`)) +
      geom_line() + xlab("Target FDR") + ylab("FDR") +ggtitle("FDR") + 
      theme(plot.title = element_text(size=22,hjust = 0.5),legend.title = element_blank()) 
    
  })
}

shinyApp(ui = ui, server = server)

