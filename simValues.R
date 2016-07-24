# This file is set up to run a loop to generate samples for 
# different sequence lengths for 1000 trials.  

source("/Users/ataylor/research/colaboration/Tasmania/R-code/treeMeasures-Unbiased.R")

star <- c(0,0, .01, .01, .01, .01, .3, .3, .3,.3)
ex1 <- c(.1, .1, .1, .1, .1, .1, .1, .1, .1, .1); 
ex2 <- c(.2, .2, .2, .2, .2, .2, .2, .2, .2, .2);
ex3 <- c(.3, .3, .3, .3, .3, .3, .3, .3, .3, .3);
Fels <- c(.05, .05, .05, .05, .3, .3, .05, .05, .3, .3);
Farris <- c(.05, .05, .05, .05, .05, .05, .3, .3, .3, .3);
bias <- c(0, 0, .05, .05, .05, .05, .2, .2, .2, .2);
bias2 <- c(0, 0, .1, .1, .1, .1, .1, .1, .1, .1); 
### The root probability.  
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

### Currently there is a bug and it is in both measures files. 
step <- 5
simLoop(Fels, .5, 5) -> data

### Graphing the results.  Noet that one needs to change only the first entry in 
### simLoop function above to use the different sets of parameters.  

title <- "Variable Sequence Length - Star Tree"
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
text(c(60,60,60),c(175, 100, 25),
     labels = c("Green: Unbiased Signed Squangles",
                "Red: Dashed Signed Minors/Solid Unsigned Minors",
                "Blue: Dashed Signed Squangle/Solid Unsigned Squangle"))

par(new=F)
