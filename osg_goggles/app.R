library(shiny)
library(tidyverse)
library(ggplot2)
library(shinyBS)

load('TidyGearCode20.Rdata')

#### Z scored abundance ####

#grab basic haul data and counts
CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, systemZone, Scientificname, N2))

#full site by species matrix
SiteXSpeciesFull <- CleanHauls %>%
  pivot_wider(id_cols = Reference:systemZone,
              names_from = Scientificname,
              values_from = N2,
              values_fill = 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  ungroup()

#collapse full site x species matrix to long term means
LTxSpeciesMeans <- SiteXSpeciesFull %>%
  subset(select = -c(Reference, systemZone, seasonYear)) %>%
  group_by(system, season) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

#make matching standard deviation matrix to long term means
LTxSpeciesStdev <- SiteXSpeciesFull %>%
  subset(select = -c(Reference, systemZone, seasonYear)) %>%
  group_by(system, season) %>%
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

#Z score conversion process
YearXSpeciesZ <- SiteXSpeciesFull %>%
  subset(select = -c(Reference, systemZone)) %>%
  group_by(system, season, seasonYear) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>% #collapse to annual means 
  ungroup() %>%
  pivot_longer(cols = !c(system:seasonYear),
               names_to = "Scientificname",
               values_to = "avg") %>%
  #expand back out to long form for leftjoins with LT mean and LT stdev
  ungroup() %>%
  left_join(pivot_longer(data = LTxSpeciesMeans, #LT mean column added
                         cols = !c(system:season),
                         names_to = "Scientificname",
                         values_to = "LTmean")) %>%
  left_join(pivot_longer(data = LTxSpeciesStdev, #LT stdev column added
                         cols = !c(system:season),
                         names_to = "Scientificname",
                         values_to = "stdev")) %>%
  group_by(system, season, seasonYear) %>%
  mutate(zscore = ((avg - LTmean)/stdev)) %>% #calculate zscores using annual means
  ungroup() %>%
  filter(LTmean > 0) %>% #remove taxa entirely absent from each system - if LTmean = 0, taxa never observed
  filter(Scientificname != "No fish") %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))

taxaList = unique(YearXSpeciesZ$Scientificname) %>%
  sort()

#### z scored hydro data #####
cleanHydro <- SiteXSpeciesFull %>%
  select(Reference:seasonYear) %>%
  left_join(HydroList)

#collapse full site x hydro matrix to long term means
LTxHydroMeans <- cleanHydro %>%
  subset(select = -c(Reference, seasonYear, Sampling_Date, year, month, Depth)) %>%
  group_by(system, season) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

LTxHydroStdev <- cleanHydro %>%
  subset(select = -c(Reference, seasonYear, Sampling_Date, year, month, Depth)) %>%
  group_by(system, season) %>%
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

#Z score conversion process
YearXHydroZ <- cleanHydro %>%
  subset(select = -c(Reference, Sampling_Date, year, month, Depth)) %>%
  group_by(system, season, seasonYear) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>% #collapse to annual means 
  ungroup() %>%
  pivot_longer(cols = !c(system:seasonYear),
               names_to = "param",
               values_to = "avg") %>%
  #expand back out to long form for leftjoins with LT mean and LT stdev
  ungroup() %>%
  left_join(pivot_longer(data = LTxHydroMeans, #LT mean column added
                         cols = !c(system:season),
                         names_to = "param",
                         values_to = "LTmean")) %>%
  left_join(pivot_longer(data = LTxHydroStdev, #LT stdev column added
                         cols = !c(system:season),
                         names_to = "param",
                         values_to = "stdev")) %>%
  group_by(system, season, seasonYear) %>%
  mutate(zscore = ((avg - LTmean)/stdev)) %>% #calculate zscores using annual means
  ungroup() %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))

hydroList = unique(YearXHydroZ$param) %>%
  sort()

######### UI ##########
# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("FWRI FIM Data Exploration"),

    # Sidebar with a slider input for number of bins 
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
        # Show a plot of the generated distribution
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
           # Show a plot of the generated distribution
           plotOutput("hydroPlot")
         )
)
))

######## server ########
# Define server logic required to draw a histogram
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
