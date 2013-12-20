library(shiny)

# Define UI
shinyUI(pageWithSidebar(
  headerPanel("Basic ODE Models of Disease Dynamics"),
  sidebarPanel(
    
    wellPanel(tags$b("Basics:"),
              sliderInput("b", "Probability of transmission:", 
                          value=.75, min=0, max=1, step=.01),
              sliderInput("k", "Average contacts (per week):", 
                          value=12, min=1, max=100, step=1)),
    
    wellPanel(tags$b("Disease properties:"),
              checkboxInput("norecov", "No recovery", value=FALSE),
              conditionalPanel(condition="input.norecov == false", 
                               sliderInput("D", "Duration (days):", 
                                           value=7, min=1, max=60, step=1),
                               sliderInput("L", label="Latent period (days):", 
                                           value=0, min=0, max=60, step=1),
                               sliderInput("seasonal", label="Seasonal fluctuations:",
                                           min=0, max=.4, step=.005, value=0))),
    
    wellPanel(tags$b("Vital dynamics and vaccination:"),
              sliderInput("bdrate", label="Birth and death rate:",
                          min=0, max=.04, step=.001, value=0),
              conditionalPanel(
                condition="input.bdrate>0",
                sliderInput("vacc", label="Proportion vaccinated at birth:",
                            min=0, max=1, step=.01, value=0),
                conditionalPanel(
                  condition="input.vacc>0",
                  sliderInput("vacceff", label="Vaccine effectiveness:",
                              min=0, max=1, step=.01, value=1)))),    
    
    wellPanel(tags$b("Time scale:"),
              selectInput("timex", label="", 
                          choices=list("Days"="days", "Years"="years")),
              conditionalPanel(condition="input.timex == 'days'", 
                               sliderInput("tmaxday", "Time max:", 
                                           value=65, min=5, max=180, step=5)),
              conditionalPanel(condition="input.timex == 'years'", 
                               sliderInput("tmaxyear", "Time max:", 
                                           value=50, min=1, max=100, step=.5))),
    
    wellPanel(tags$b("Plot settings:"),
              #                   checkboxInput("maxinfect", "Hide infection info", value=FALSE),
              checkboxInput("inf.only", "Show infectious curve only", value=FALSE),
              conditionalPanel(
                condition="input.bdrate>0",
                checkboxInput("second.inf", "Post-initial outbreak only", value=FALSE))),
    
    wellPanel(tags$h5("Created by Mathew Kiang"), 
              tags$body("(", tags$a("Git", 
                                    href="http://github.com/mkiang")," | ", 
                        tags$a("Tweet", 
                               href="http://twitter.com/mattkiang")," | ",
                        tags$a("Read", 
                               href="http://mathewkiang.com"),")"))),
  
  mainPanel(plotOutput("guessPlot"),
            wellPanel(
              tags$body(h3("About"), 
                        p(("This is just a tool to easily visualize some basic ODE models—specifically, SI, SIR, and SIER models. You can add vital dynamics (however, birth and death rates are equal to maintain a constant population) and vaccination programs (incorportaing vaccine effectiveness). It's certainly not a comprehensive set of all math models of disease. I'm hoping to find time to code an analogous (stochastic) network-based disease dynamics page."), 
                          h3("Things to note"), 
                          p("Not all permutations of every parameter will render a plot (or a sensible plot). This is most often a function of a time scale that doesn't make sense given the parameters (e.g., taking a basic SIR model out to 100 years), so try adjusting your time scale first. I don't perform any computational adjustments to the time scale so increasing your time frame will result in (potentially) exponentially increased waiting times. The model always assumes a constant population of 1,000,000 people and always starts with a single infected person in an otherwise naive population. Lastly, vaccinations only occur at birth."),
                          p("The code could certainly be extended to incorporate more complex features (e.g., disease vectors, age structures, etc.), and I encourage you to use my code on Github to do so."))
              ))
  )))
