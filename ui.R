library(shiny)

# Define UI
shinyUI(pageWithSidebar(
  headerPanel("Disease Dynamics using ODE Models"),
  sidebarPanel(
    
    wellPanel(tags$b("Basics:"),
      sliderInput("b", "Probability of transmission:", 
                  value=.75, min=0, max=1, step=.01),
      sliderInput("k", "Average contacts (per week):", 
                  value=12, min=1, max=100, step=1),
      sliderInput("D", "Disease duration (days):", 
                  value=7, min=1, max=60, step=1),
      sliderInput("L", label="Incubation period (days):", 
                  value=0, min=0, max=60, step=1)),
    
    wellPanel(tags$b("Time control:"),
      selectInput("timex", label="", 
                  choices=list("Days"="days", "Months"="months")),
      conditionalPanel(condition="input.timex == 'days'", 
                       sliderInput("tmaxday", "Time max:", 
                       value=65, min=5, max=365, step=5)),
      conditionalPanel(condition="input.timex == 'months'", 
                       sliderInput("tmaxmonth", "Time max:", 
                       value=600, min=6, max=900, step=6))),
    
    wellPanel(tags$b("Advanced:"),
      sliderInput("bdrate", label="Birth and death rate:",
                  min=0, max=.03, step=.001, value=0),
      sliderInput("seasonal", label="Seasonal fluctuations:",
                  min=0, max=.3, step=.005, value=0),
      checkboxInput("maxinfect", "Hide infection info", value=FALSE),
      checkboxInput("logy", "Use log for y-axis", value=FALSE),
      checkboxInput("inf.only", "Show infections only", value=FALSE)),
    wellPanel(tags$h5("Created by Mathew Kiang"), 
              tags$body("(", tags$a("Git", 
                                    href="http://github.com/mkiang")," | ", 
                        tags$a("Tweet", 
                               href="http://twitter.com/mattkiang")," | ",
                        tags$a("Read", 
                               href="http://mathewkiang.com"),")"))),
  
  mainPanel(plotOutput("guessPlot"),
            wellPanel(
              tags$body(h3("About the models"), 
                        p(("These are simple Susceptible-Infected-Recovered (SIR) or Susceptible-Exposed-Infected-Recovered (SEIR) models that can be used to simulate an epidemic within a large population. For more information about compartmental models in general, see"), a("the Wikipedia page.", href="http://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology")), p("For more information about these models specifically, including assumptions, equations, and a more thorough explanation of the parameters, see my original post (INSERT LINK)."), h3("Things to note"), p("For simplicity, the model always assumes the epidemic starts with just one infectious person in a constant population of 1 million people with birth rate equal to death rate. By setting the birth and death rate to zero, you effectively have a closed population with drastically different diseaase dynamics than an open population. Switching to an SIR from an SEIR is as simple as setting the incubation period to 0 (the legend at the bottom of the graph will help you keep track of which model you are using). Switch the time scale to months to see long term population dynamics (seasonal effects in an SEIR model with an open population need much longer time scales to reveal themselves.)"),
                        p("About performance: To ensure reponsiveness, the larger time max becomes, the less fine-grained (i.e., less accurate) it is and the longer computation will take. Further, each ticking or unticking prompts the server to rerun the entire model--you'll need to be patient or lower the value of time max."), 
                        p("About accuracy: I'm fairly confident most permutations of these variables will be relatively accurate. That said, the models were modified from my homework assignments--sometimes very heavily modified--to accommodate compressing many different models into one workable web-app. They were not doublechecked. There is definitely something going on with the seasonal effects over very long time periods, but I haven't dedicated the time to figuring out what that is and I probably won't since it is accurate enough to give you an idea of what is supposed to happen."),
                        p("Enjoy."))
            ))
))
