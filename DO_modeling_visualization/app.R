library(shiny)
library(shinythemes)
library(tidyverse)
library(readr)
library(lubridate)
start_stop <- read.csv("start_stop.csv")
start = as.Date(start_stop$start_all)
stop = as.Date(start_stop$stop_all)
year<-2019

# Define UI
ui <- fluidPage(theme = shinytheme("lumen"),
                titlePanel(h1("Dissolved Oxygen Forecast",align = "center")),
                title = "DO forecast",
                sidebarPanel(
                  selectInput(inputId = "model", label = strong("Model Type"),
                              choices = c("SSS only" ="_sss","All" = "_all","Oxygen only" = "_o2","Temperature only" = "_temp"),
                              selected = "_all"),
                  selectInput(inputId = "year", label = strong("Forecast Year"),
                              choices = c("2019","2018","2017","2016","2015","2014","2013"),
                              selected = "2019"),
                  
                  conditionalPanel(
                    condition = "input.gif == false",
                    sliderInput("date", strong("Forecast date"), value = min(start[year(start)==year]),
                                min = min(start[year(start)==year]), max = max(stop[year(stop)==year]), step = 1)
                  ),
                  checkboxInput("param","Display parameter evolution", value = F),
                  checkboxInput("gif","Display forecasts as gif", value = F),
                  checkboxInput("eval","Display forecast evaluation", value = F)
                ),
                
                # Output: 
                mainPanel(
                  imageOutput("preImage"),br(),br(),
                  textOutput("forecast_date"),br(),
                  conditionalPanel(
                    condition = "input.param == true",
                    imageOutput("param")
                  ),
                  conditionalPanel(
                    condition = "input.eval == true",
                    h3("Forecast evaluation: RMSE"),
                    imageOutput("rmse"),
                    h3("Forecast evaluation: Uncertainty Analysis"),
                    h4("Partitioned uncertainty over time"),
                    imageOutput("uncert_year"),br(),br(),
                    h4("Partitioned uncertainty as a function of forecast horizon"),
                    imageOutput("uncert_all")
                  )
                )
)

#Define server function
server <- function(input, output, session) {
  # Send a pre-rendered image, and don't delete the image after sending it
  output$preImage <- renderImage({
    year <- input$year
    model <- input$model
    date <- as.Date(input$date)
    start_thisYear = as.Date(start_stop$start_all[year(start_stop$start_all) == year])
    date_num = formatC(as.numeric(difftime(date,start_thisYear))+1, width = 3, format = "d", flag = "0")
    if(model == "_all"){model <- ""}
    filename <- paste('06May20/',year,model,'/forecast', date_num, '.png', sep='')
    if(input$gif == T){
      filename <- paste('06May20/',year,model,'/animated_forecast.gif',sep="")
    }
    # Return a list containing the filename and alt text
    list(src = filename,
         alt = paste("Forecast failed to load"))
  }, deleteFile = FALSE
  )
  
  
  output$forecast_date <- renderText({
    paste("Forecast generated on 06 May 2020")
  })
  
  
  output$rmse <- renderImage({
    filename <- "./06May20/RMSE_all.jpeg"
    height <- session$clientData$output_preImage_height
    list(src = filename,
         height = height,
         alt = paste("RMSE failed to load"))
  }, deleteFile = F)
  
  
  output$uncert_all <- renderImage({
    filename <- "./06May20/uncert_horiz.png"
    height <- session$clientData$output_preImage_height
    list(src = filename,
         height=height,
         alt = paste("Uncertainty figure failed to load"))
  }, deleteFile = F)
  
  
  output$uncert_year <- renderImage({
    filename <- "./06May20/uncert_season7.png"
    height <- session$clientData$output_preImage_height
    list(src = filename,
         height=height,
         alt = paste("Second uncertainty figure failed to load"))
  }, deleteFile = F)
  
  output$param <- renderImage({
    model <- input$model
    if(model == "_all"){model <- ""}
    filename <- paste('06May20/',input$year,model,'/params.jpeg', sep='')
    height <- session$clientData$output_preImage_height
    list(src = filename,
         height=height,
         alt = paste("Parameter evolution figure failed to load"))
  }, deleteFile = F)
  
  output$forecast_date <- renderText({
    paste("Forecast generated on 06 May 2020")
  })
  
  
  observe({
    updateSliderInput(session, "date", value = min(start[year(start)==input$year]), min = min(start[year(start)==input$year]), 
                      max = max(stop[year(stop)==input$year]))
  })
}

# Create Shiny object
shinyApp(ui = ui, server = server)