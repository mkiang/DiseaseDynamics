DiseaseDynamics
===============

## Introduction

This is my current project to use R, Shiny, and deSolve and create simple ODE models for disease dynamics (e.g., SIR, SEIR). 
The other projects (BasicSIR, BasicSEIR, OpenSEIR) were all starting points to get familiar with Shiny, but this project should replace them in both functionality and priority.
To see the code in action, go to: See the code in action at: https://mkiang.shinyapps.io/DiseaseDynamics/

Short blog post describing the parameters: http://mathewkiang.com/2013/12/20/shiny-desolve-interactive-ode-models/


### To run the model on your own computer
If you want to run this on your own computer, make sure you have `shiny`, `quantmod`, and `deSolve` installed. Then run the following command: `shiny::runGitHub('DiseaseDynamics', 'mkiang')`

Or copy and paste the following script.

```r
## Running shiny() on your own computer
## Install libraries -- shiny, quantmod, deSolve
if (!require(shiny)) {
    install.packages("shiny")
    library(shiny)
}
if (!require(quantmod)) {
    install.packages("quantmod")
    library(quantmod)
} 
if (!require(deSolve)) {
    install.packages("deSolve")
    library(deSolve)
} 

## Load libraries
require(deSolve)
require(quantmod)
require(shiny)

## Run shinyapp
shiny::runGitHub('DiseaseDynamics', 'mkiang')
```
