# All of these models are (loosely) derived from assignments in my Dynamics of 
# Infectious Disease (EPI 501) course at HSPH. There is a nonzero probability
# errors are contained within. 

## Load packages we'll need.
library(shiny)
library(deSolve) 
library(ggplot2)

## Define various comparmental models we will use later. All models are in an
## open population, but that changes by simply setting births and deaths to 0.
#   1) N is total population. S is susceptible. E is exposed (but not 
#       infectious). I is infectious. R is recovered. 
#   2) b is the probability of transmission from one infected to a susceptible
#       at one point in time (i.e., one contact).
#   3) k is the contact rate.
#   4) Models that incorporate seasonal variation decompose k such that:
#       k = k0 * (1 + k1 * cos(2 * pi * t)) 
#         where k0 is the overall average contact rate, k1 is the amplitude
#         (fluctuation about k0) and t is time. For simplicity, this can only
#         be added when the time scale is yearly.
#   5) Together, bk = beta in our models.
#   6) r is the rate of recovery--also the reciprocal of disease duration.
#   7) L is the latent period where a person is neither infectious nor 
#       susceptible.

### Basic SIR open population 
SIR.open <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    N <- S + I + R
    dS <- birth*N - beta*S*I/N - death*S
    dI <- beta*S*I/N - r*I - death*I
    dR <- r*I - death*R
    result <- c(dS, dI, dR)
    list(result)
  })
}

### Basic SEIR open population
SEIR.open <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    N <- S + E + I + R
    dS <- birth*N - beta*S*I/N - death*S
    dE <- beta*S*I/N - L*E - death*E
    dI <- L*E - r*I - death*I
    dR <- r*I - death*R
    result <- c(dS, dE, dI, dR)
    list(result)
  })
}

### SEIR with seasonal variation
SIR.kt <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    k <- k0 * (1 + k1 * cos(2 * pi * (t/t.adj)))  
    beta <- b * k
    N <- S + I + R
    dS <- -beta*S*I/N + birth*N - death*S
    dI <- beta*S*I/N - r*I - death*I
    dR <- r*I - death*R
    result <- c(dS, dI, dR)
    list(result)
  })
}

### SEIR with seasonal variation
SEIR.kt <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    k <- k0 * (1 + k1 * cos(2 * pi * (t/t.adj)))  
    beta <- b * k
    N <- S + E + I + R
    dS <- -beta*S*I/N + birth*N - death*S
    dE <- beta*S*I/N - k*E - death*E
    dI <- k*E - r*I - death*I
    dR <- r*I - death*R
    result <- c(dS, dE, dI, dR)
    list(result)
  })
}

## Define a function to find max infections to put on plots later.
max.inf <- function(x) {
  infectious <- round(max(x$people[x$compartments=="I"], na.rm=TRUE))
  timestep <- which.max(x$people[x$compartments=="I"])
  info <- cbind(infectious, timestep)
  return(info)
}

## Makes the data in a more ggplot-friendly format
stacker <- function(df){
  df.stack <- stack(df[, -1])
  df.stack$time <- rep(seq_len(nrow(df)), length(table(df.stack$ind)))
  names(df.stack)[1:2] <- c("people", "compartments")
  df.stack$compartments <- factor(df.stack$compartments, 
                                  levels=c("S","E","I","R"), ordered=TRUE)
  return(df.stack)
}

## Start the Shiny stuff
shinyServer(function(input, output) { 
  output$guessPlot <- renderPlot(function() {
    
    # To optimize responsiveness (i.e., minimize computation), we change
    # the scale when the number of time points gets too high. However,
    # this means we also need to change the x-axis labels and the annotation
    # on the graph. We will use divtime to do that. We also need to know
    # if we are using months or days in order to adjust the rates.
    if (input$timex=="days") {
      r.adj <- 1; k.adj <- 1/7; bd.adj <- 365; L.adj <- 1; t.adj <- 365
      timemax <- input$tmaxday
    } else {
      r.adj<-365/12; k.adj<-52/12; bd.adj<-12; L.adj<-365/12; t.adj<-12
      timemax <- input$tmaxmonth
    }
    if (timemax <= 180) divtime <- 100 else divtime <- 10
    dt <- seq(0, timemax, by=1/divtime)
    
    # First, making sure we use the right model.
    ## SIR without seasonal effects
    if(input$L==0 & input$seasonal==0)  {
      r <- (1/input$D) * r.adj
      b <- input$b
      k <- input$k * k.adj
      birth <- death <- input$bdrate/bd.adj
      
      inits <- c(S=1e6-1, I=1, R=0)
      guess_pars <- c(r=r, beta=b*k, birth=birth, death=death)
      guess <- stacker(as.data.frame(ode(inits, dt, SIR.open, 
                                         parms=guess_pars)))
    } else if (input$L!=0 & input$seasonal==0) {
   
    ## SEIR without seasonal effects
      r <- (1/input$D) * r.adj
      b <- input$b
      k <- input$k * k.adj
      L <- (1/input$L) * L.adj
      birth <- death <- input$bdrate/bd.adj
      
      inits <- c(S=1e6-1, E=0, I=1, R=0)
      guess_pars <- c(r=r, beta=b*k, birth=birth, death=death, L=L)
      guess <- stacker(as.data.frame(ode(inits, dt, SEIR.open, 
                                         parms=guess_pars)))
    } else if (input$L==0 & input$seasonal!=0)  {
    
    ## SIR with seasonal effects
      r <- (1/input$D) * r.adj
      b <- input$b
      k0 <- input$k * k.adj
      k1 <- input$seasonal
      birth <- death <- input$bdrate/bd.adj
      
      inits <- c(S=1e6-1, I=1, R=0)
      guess_pars <- c(r=r, b=b, birth=birth, death=death, k0=k0, k1=k1,
                      t.adj=t.adj)
      guess <- stacker(as.data.frame(ode(inits, dt, SIR.kt, 
                                         parms=guess_pars)))
    } else if (input$L!=0 & input$seasonal!=0)  {
      
    ## SEIR with seasonal effects
      r <- (1/input$D) * r.adj
      b <- input$b
      k0 <- input$k * k.adj
      k1 <- input$seasonal
      L <- (1/input$L) * L.adj
      birth <- death <- input$bdrate/bd.adj
      
      inits <- c(S=1e6-1, E=0, I=1, R=0)
      guess_pars <- c(r=r, b=b, birth=birth, death=death, k0=k0, k1=k1,
                      t.adj=t.adj, L=L)
      guess <- stacker(as.data.frame(ode(inits, dt, SEIR.kt, 
                                         parms=guess_pars)))
    }
    
    # Max infections
    max.guess <- as.data.frame(max.inf(guess))
    max.guess$compartments <- "I"
    max.guess$adj.timestep <- max.guess$timestep / divtime
    
    # Base plots
    if (input$inf.only==FALSE) {
    base.plot <- ggplot(guess, aes(x=time, y=people, 
                                   group=compartments, color=compartments)) + 
      xlab("Time") + ylab("Number of People")
    } else {
    base.plot <- ggplot(guess[guess$compartments=="I",], aes(x=time, y=people, 
                                   group=compartments, color=compartments)) + 
      xlab("Time") + ylab("Number of Infected")
    }
    label.plot <- base.plot + geom_line(alpha=.8) + 
      theme(legend.position="bottom") + theme(legend.title=element_blank())
    
    # Conditioning plots
    if (input$maxinfect==TRUE) {
      rescaledx.plot <- label.plot + 
        scale_x_continuous(labels=function(x) x/divtime) 
    } else {
        rescaledx.plot <- label.plot + 
        scale_x_continuous(labels=function(x) x/divtime) + 
        geom_text(data=max.guess, 
                  aes(timestep*1.2, infectious, 
                      label=paste(infectious, "infections\nat", adj.timestep)), 
                  hjust=0, vjust=1, show_guide=FALSE, size=5, 
                  alpha=.5, color="black")
    }
    
    # Log(y)?
    if (input$logy==TRUE) {
      final.plot <- rescaledx.plot + coord_trans(y="log10")
    } else {
      final.plot <- rescaledx.plot
    }
    
    
    print(final.plot)
    
  })  # closes reactivePlot
})  # Closes shinyServer
