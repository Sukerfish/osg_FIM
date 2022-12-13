#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(ggplot2)

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

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("FWRI FIM Data"),

    # Sidebar with a slider input for number of bins 
    wellPanel(
          selectInput(
            "display_var",
            "Which taxa to display?",
            choices = c(taxaList)
          
        )),
wellPanel(
        # Show a plot of the generated distribution
           plotOutput("taxaPlot")
        )
    )


# Define server logic required to draw a histogram
server <- function(input, output) {

  chunk <- eventReactive(input$display_var, {
    filter(YearXSpeciesZ, Scientificname %in% c(input$display_var))
  })
  
  #science plot
  output$taxaPlot <- renderPlot({
    ggplot(data = chunk(),#dynamically filter the sci variable of interest
           aes(x=seasonYear,
               y=zscore)) +
      geom_point() +
      geom_line() +
      theme(axis.text=element_text(size = 12)) +
      theme(axis.title=element_text(size = 16)) +
      theme(strip.text = element_text(size = 16)) +
      facet_grid(season~system)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
