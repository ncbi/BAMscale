list.of.packages <- c("shiny", "ggplot2", "tidyr", "ggrepel", "gridExtra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
                      if(length(new.packages)) install.packages(new.packages)
                      

library(shiny)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(gridExtra)

options(shiny.maxRequestSize=150*1024^2)  # Current Max File Size: 150MB  (For larger files increase this accordingly)

# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  ),
  fluidRow(
    column(
      6, 
      align = "center",
      tags$img(src = "https://healthtech.upenn.edu/wp-content/uploads/2018/09/nci_-1200x600.jpg", align = "left", height = "100px", width = "200px")   #NIH-NCI
    )
  ),
  fluidRow(
    column(12,
           align = "center",
           titlePanel(
             title  = "Application for Plotting Peak-Density Figures",
             windowTitle = "Application for Plotting Peak Density Figures"
             ),
           fluidRow(
             column(
               3,
               align = "center",
               fileInput(
                 "file1", 
                 "Choose Normalized Coverage File",
                 accept = c(".tsv")),
               tags$hr()
             ),
             column(
               6,
               align = "center",
               uiOutput("slider")
             )
           ),
           fluidRow(
             column(
               6,
               align = "center",
               uiOutput("xAxis")
             ),
             column(
               6, 
               align = "center",
               uiOutput("yAxis")
             )
           ),
           fluidRow(
             column(
               12,
               align = "center",
               plotOutput("plotXY", width = "700px", height = "700px")
             )
           ),
           fluidRow(
             column(
               12,
               align = "right",
               "Developed by Jacob M Gross & Lorinc S Pongor"
             )
           )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  raw_data = eventReactive(input$file1, {
    a = read.table(input$file1$datapath, stringsAsFactors = F, check.names = F, sep = "\t", header = T, fill = T)
    a = separate(a, 1, c("chr", "position"), sep = ":")
    a = separate(a, 2, c("Start", "Stop"), sep = "-")
  })
  
  output$slider = renderUI({
    dat = as.data.frame(t(rbind(raw_data()[, input$xaxis], raw_data()[, input$yaxis])))
    sliderInput(
      "slider1", 
      label = h3("Axis Limit Slider"), 
      min = 0, 
      max = round(max(dat)+1), 
      value = max(quantile(dat[,1], .99), quantile(dat[,2], .99)) 
    )
  })
  
  output$xAxis = renderUI({
    dat1 = as.data.frame(raw_data()[,4:ncol(raw_data())])
    items = names(dat1)
    selectInput("xaxis", "X-axis:", items, selected = items[1])
  })
  
  output$yAxis = renderUI({
    dat1 = as.data.frame(raw_data()[,4:ncol(raw_data())])
    items = names(dat1)
    #names(items) = items
    selectInput("yaxis", "Y-axis:", items, selected = items[2])
  })
   
   output$plotXY <- renderPlot({
     dat = as.data.frame(t(rbind(raw_data()[, input$xaxis], raw_data()[, input$yaxis])))
     colnames(dat) = c(paste0("1: ", input$xaxis), paste0("2: ", input$yaxis))
     
     colpal = rev(rainbow(5))
     colpal[1] = "#562188FF"
     colpal[3] = "#225A16FF"
     
     p1 = ggplot(dat, aes(x = dat[,1], y = dat[,2])) +
       labs(x = paste0(colnames(dat)[1], " (Reads Per Peak)") , y = paste0(colnames(dat)[2], " (Reads Per Peak)")) +
       ggtitle(label = paste0(colnames(dat)[1]," vs ",colnames(dat)[2], " LibSize Filtered")) +
       geom_hex(bins = 200) +
       scale_fill_gradientn(colours = colpal, name = "log2\n(Density of Peaks)", trans = "log2") +
       scale_x_continuous(expand = c(0, 0), limits = c(0, input$slider1)) +
       scale_y_continuous(expand = c(0, 0), limits = c(0, input$slider1)) +
       geom_abline(intercept = 0, slope = c(0.15579, 0.3193, 0.7208, 1.38742, 3.1315, 6.4188), colour="#8B8B8B", size = 0.5) +
       geom_abline(intercept = 0, slope = c( 0.5, 2), colour="#7E7D7D", size = 0.85)+
       geom_abline(intercept = 0, slope = 1, colour="black", size=0.8, linetype = "twodash") +
       scale_colour_gradient(low = "red", high = "blue") +
       coord_equal(ratio = 1) +
       theme(plot.title = element_text(size = 10), legend.title = element_text(size=8), axis.title = element_text(size=8.5),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "#F3F3F3", colour = "black"), 
             axis.text = element_text(size = 16), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))
     
     grid.arrange(p1)
   })
 
}

# Run the application 
shinyApp(ui = ui, server = server)

