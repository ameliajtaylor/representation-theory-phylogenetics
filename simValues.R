# This file is set up to run a loop to generate samples for 
# different sequence lengths for 1000 trials.  
# The parameters below match the ones Barbara ran for the paper. 

source("/Users/ataylor/research/colaboration/Tasmania/R-code/treeMeasures-Unbiased.R")

star <- c(0,0, .01, .01, .01, .01, .3, .3, .3,.3)
ex1 <- c(.1, .1, .1, .1, .1, .1, .1, .1, .1, .1); 
ex2 <- c(.2, .2, .2, .2, .2, .2, .2, .2, .2, .2);
ex3 <- c(.3, .3, .3, .3, .3, .3, .3, .3, .3, .3);
Fels <- c(.05, .05, .05, .05, .3, .3, .05, .05, .3, .3);
Farris <- c(.05, .05, .05, .05, .05, .05, .3, .3, .3, .3);
bias <- c(0, 0, .05, .05, .05, .05, .2, .2, .2, .2);
bias2 <- c(0, 0, .1, .1, .1, .1, .1, .1, .1, .1); 
### The root distribution. 
p <- .4
###sequence length is N. Trials is R. 

### The loop runs 
simLoop <- function(tree, p, step){
  resultArray <- NULL
  for (i in seq(10, 100, step)){
    sample <- F12.34(tree, p, i, 1000, 2)
    SSSm <- measuresRSS(sample, i, 1000, stRSS) 
    USSm <- measures(sample, st.measure, i, 1000)
    SSSs <- measuresRSS(sample, i, 1000, markovRSS) 
    USSs <- measures(sample, markovInv, i, 1000)
    UnSSS <- measuresUnBiased(sample, i, 1000) 
    resultArray <- c(resultArray, c(SSSm[1], USSm[1], SSSs[1], USSs[1], UnSSS[1]))
  }
  matrix(resultArray, 5, 100/step - 1)
}

step <- 5
simLoop(Fels, .5, 5) -> data
### Below we build plots of all the different measures.  Note that to run different 
### parameters, all we need to change is the first argument of the simLoop function above. 

### Currently there is a bug.  I thought it was in the unbiased squangle, but now 
### believe it is elsewhere. 

title <- "Variable Sequence Length - Felsenstein Zone"
  plot(seq(10,100,step),data[1,], main=title, 
         xlab="sequence length", ylab="number of success in 1000",
       pch =".", xlim=c(0,100), ylim=c(0,1000), col="red")
  par(new=T)
  plot(seq(10,100,step),data[2,], main=title, 
      xlab="sequence length", ylab="number of success in 1000",  
       pch =".", xlim=c(0,100), ylim=c(0,1000), col="red")
  par(new=T)
  plot(seq(10,100,step),data[3,], main=title, 
      xlab="sequence length", ylab="number of success in 1000",  
       pch =".", xlim=c(0,100), ylim=c(0,1000), col="blue")
  par(new=T)
  plot(seq(10,100,step),data[4,], main=title, 
       xlab="sequence length", ylab="number of success in 1000",  
       pch =".", xlim=c(0,100), ylim=c(0,1000), col="blue")
  par(new=T)
  plot(seq(10,100,step),data[5,], main=title, 
     xlab="sequence length", ylab="number of success in 1000",  
     pch =".", xlim=c(0,100), ylim=c(0,1000), col="green")
  lines(seq(10,100,step), data[1,], col="red", lty=2) #SSSminors
  lines(seq(10,100,step), data[2,], col="red", lty=1) #USSminors
  lines(seq(10,100,step), data[3,], col="blue", lty=2) #SSSsquangle
  lines(seq(10,100,step), data[4,], col="blue", lty=1) #USSsquangle
  lines(seq(10,100,step), data[5,], col="green",lty=2) #Unbiased Signed Squangle
text(c(100),c(250), 
     labels = ('Green: Unbiased Squangle'))

  par(new=F)
  
  c(1000, 925, 850, 775)
  # 100, 100, 100    , 175, 100, 25
  #"Dashed Red Signed Minors", 
  #"Solid Blue Unsigned Squangle", "Dashed Blue Signed Squangle"