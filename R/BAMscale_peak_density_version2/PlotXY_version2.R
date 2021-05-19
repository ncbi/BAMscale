library(shiny)
library(GenomicRanges)
library(data.table)
library(ggplot2)

options(shiny.maxRequestSize=300*1024^2)


subsetQuantwithBed = function(peaksfile, bedfile) {
  coords = read.table(bedfile)[,1:3]
  colnames(coords) = c("chr", "start", "end")
  coords = makeGRangesFromDataFrame(coords)
  
  peaks = read.table(peaksfile,
                     header = T,
                     sep = "\t")
  
  peakcoords = GRanges(peaks[,1])
  
  overs = findOverlaps(query = peakcoords, subject = coords)
  peakoverlaps = unique(queryHits(overs))
  
  return (peaks[peakoverlaps,])
  
}


ui <- fluidPage(
  fluidRow(
    column(
      6, 
      align = "center",
      tags$img(src = "https://healthtech.upenn.edu/wp-content/uploads/2018/09/nci_-1200x600.jpg", align = "left", height = "100px", width = "200px")   #NIH-NCI
    )
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput('datafile', 'Choose quantified peak file',
                accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
      
      fileInput("bedfile", "Choose a BED file to subset peaks",
                multiple = FALSE,
                accept = c(".bed")),

      conditionalPanel(
        # use a server side condition
        # placeholders will be replaced from the server
        condition = "output.fileUploaded",
        selectInput("xAxis", "X-axis sample", ""),
        selectInput("yAxis", "Y-axis sample", ""),
        numericInput("slidelimit", "Axis limit", 1),
        numericInput("hexcount", "No. of hex bins", 200))
      
      
      
    ),
    mainPanel(
      plotOutput("plotXY")
    )
  )
)

server <- function(input, output, session){
  # create reactive version of the dataset (a data.frame object)
  filedata <- reactive({
    infile <- input$datafile
    if (is.null(infile))
      # User has not uploaded a file yet. Use NULL to prevent observeEvent from triggering
      return(NULL)
    
    if(is.null(input$bedfile)) {
      temp <- read.table(input$datafile$datapath,
                       header = T,
                       sep = "\t", check.names = F)
    } else {
      temp = subsetQuantwithBed(input$datafile$datapath, input$bedfile$datapath)
    }
  })
  
  output$fileUploaded <- reactive({
    return(!is.null(filedata()))
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  
  observeEvent(filedata(), {
    snames = names(filedata())
    snames = snames[2:length(snames)]
    updateSelectInput(session, "xAxis", choices =snames, selected=snames[1])
    updateSelectInput(session, "yAxis", choices =snames, selected = snames[2])
    xv = input$yAxis
    yv = input$xAxi
    if (is.null(xv) || is.null(yv)){
      xv = snames[1]
      yv = snames[2]
    } 
    
    df = as.data.frame(t(rbind(filedata()[, xv], filedata()[, yv])))
    axislim = round(min(quantile(df[,1], .99)[1], quantile(df[,2], .99)[1]))
    axismax = max(max(df[,1]), max(df[,2]))
    baxislim =  max(quantile(df[,1], .999)[1], quantile(df[,2], .999)[1])
    updateNumericInput(session, "slidelimit", value = axislim)
  })
 
  output$plotXY <- renderPlot({
    colpal = rev(rainbow(5))
    colpal[1] = "#562188FF"
    colpal[3] = "#225A16FF"
    
    gp1 <- NULL
    gp2 = NULL
    
    dat = filedata()
    if (!is.null(dat)){
      xv <- input$xAxis
      yv <- input$yAxis
      if (!is.null(xv) & !is.null(yv)){
        df = as.data.frame(t(rbind(filedata()[, input$xAxis], filedata()[, input$yAxis])))
        colnames(df) = c(input$xAxis, input$yAxis)
        axislim = min(quantile(df[,1], .99)[1], quantile(df[,2], .99)[1])
        axismax = max(max(df[,1]), max(df[,2]))
        baxislim =  max(quantile(df[,1], .999)[1], quantile(df[,2], .999)[1])
        
        gp1 = ggplot(df, aes(x = df[,1], y = df[,2])) +
          geom_hex(bins = input$hexcount) +
          scale_fill_gradientn(colours = colpal, name = "log2\n(Density of Peaks)", trans = "log2") +
          scale_x_continuous(expand = c(0, 0), limits = c(0, input$slidelimit)) +
          scale_y_continuous(expand = c(0, 0), limits = c(0, input$slidelimit)) +
          xlab(input$xAxis) +
          ylab(input$yAxis) +
          coord_equal(ratio = 1) +
          geom_abline(intercept = 0, slope = 1, colour="black", size=1.00, linetype = "twodash") +
          theme_classic() +
          theme(text = element_text(size = 14))
        }
    }
    
    return (gp1)
    })
}

shinyApp(ui, server)
