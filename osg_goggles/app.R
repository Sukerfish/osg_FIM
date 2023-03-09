library(shiny)
library(tidyverse)
library(ggplot2)
library(shinyBS)

load('ShinyGearCode20.RData')

taxaList = unique(YearXSpeciesZ$Scientificname) %>%
  sort()

hydroList = unique(YearXHydroZ$param) %>%
  sort()

######### UI ##########
# Define UI for application
ui <- fluidPage(

    # Application title
    titlePanel("FWRI FIM Data Exploration"),

    # control buttoms 
    tabsetPanel(
      tabPanel(title = "Biology",
               actionButton(inputId = "reload", 
                            label = "Random", 
                            icon("arrows-rotate"),
                            style = "color: #fff; background-color: #963ab7; border-color: #2e6da4"),
               actionButton(inputId = "viewchunk", 
                            label = "View Data Table",
                            style = "color: #fff; background-color: #103de2; border-color: #2e6da4"),
               bsModal("modalExample", "Data Table", "viewchunk", size = "large",
                       dataTableOutput("tbl")),
                            #icon("arrows-rotate"),
                            #style = "color: #fff; background-color: #963ab7; border-color: #2e6da4"),
      selectizeInput(
            inputId = "bio_var",
            label = "Which taxa to display?",
            choices = c(taxaList),
            selected = sample(taxaList, 1),
            multiple = TRUE,
            options = list(plugins= list('remove_button'))
        ),
      checkboxInput(
        inputId = "taxOpt",
        label = "Smooth"
      ),
wellPanel(
        # Show a plot
           plotOutput("taxaPlot")
        )
    ),
tabPanel(title = "Hydrology",
         selectizeInput(
           inputId = "hydro_var",
           label = "Which hydro data to display?",
           choices = c(hydroList),
           selected = "Temperature",
           multiple = TRUE,
           options = list(plugins= list('remove_button'))
         ),
         checkboxInput(
           inputId = "hydOpt",
           label = "Smooth"
         ),
         wellPanel(
           # Show a plot
           plotOutput("hydroPlot")
         )
)
))

######## server ########
# Define server logic 
server <- function(input, output, session) {

  chunk <- eventReactive(input$bio_var, {
    filter(YearXSpeciesZ, Scientificname %in% c(input$bio_var))
  })
  
  observeEvent(input$reload, {
  updateSelectizeInput(session, "bio_var", choices = c(taxaList), selected = sample(taxaList, 1))
  })
    
 output$tbl = renderDataTable(select(chunk(), Scientificname, LTmean, stdev, system, season) %>%
                                unique(),
                              options = list(order = list(list(0, 'asc'))))
  
  #science plot
  output$taxaPlot <- renderPlot({
    p <- ggplot(data = chunk(),#dynamically filter the sci variable of interest
           aes(x=seasonYear,
               y=zscore,
               color=Scientificname)) +
      geom_point() +
      geom_line() +
      theme(axis.text=element_text(size = 12)) +
      theme(axis.title=element_text(size = 16)) +
      theme(strip.text = element_text(size = 16)) +
      facet_grid(season~system)
    
    if (input$taxOpt){
      p <- p + geom_smooth(method=lm, se=FALSE)
    }
    
    print(p)
  })
  
  hunk <- eventReactive(input$hydro_var, {
    filter(YearXHydroZ, param %in% c(input$hydro_var))
  })
  
  #hydro plot
  output$hydroPlot <- renderPlot({
    q <- ggplot(data = hunk(),#dynamically filter the sci variable of interest
           aes(x=seasonYear,
               y=zscore,
               color=param)) +
      geom_point() +
      geom_line() +
      theme(axis.text=element_text(size = 12)) +
      theme(axis.title=element_text(size = 16)) +
      theme(strip.text = element_text(size = 16)) +
      facet_grid(season~system)
    
    if (input$hydOpt){
      q <- q + geom_smooth(method=lm, se=FALSE)
    }
    
    print(q)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
