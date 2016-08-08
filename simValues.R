# This file is set up to run a loop to generate samples for 
# different sequence lengths for 1000 trials.  

source("/Users/ataylor/research/colaboration/Tasmania/R-code/treeMeasures-Unbiased.R")

# the probabilities here are internal edge a, b, M1, M2, M3, M4.  
# We consider the root at an internal node, so a single internal edge.  

star <- c(0,0, .01, .01, .01, .01, .3, .3, .3,.3)
ex1 <- c(.1, .1, .1, .1, .1, .1, .1, .1, .1, .1); 
ex2 <- c(.05, .05, .15, .15, .15, .15, .15, .15, .15, .15);
ex3 <- c(.05, .05, .3, .3, .3, .3, .3, .3, .3, .3);
Fels <- c(.05, .05, .05, .05, .3, .3, .05, .05, .3, .3); #matches Barbara
Farris <- c(.05, .05, .05, .05, .05, .05, .3, .3, .3, .3); #matches Barbara
bias <- c(0, 0, .05, .05, .05, .05, .2, .2, .2, .2); #She does "star Farris" below.
StarFarris <- c(0, 0, .05, .05, .05, .05, .3, .3, .3, .3); 
### The root probability.  
p <- .4
###sequence length is N. Trials is R. 

### The loop runs 
simLoop <- function(tree, p, step,seqlen){
  resultArray <- NULL
  for (i in seq(10, seqlen, step)){
    sample <- F12.34(tree, p, i, 1000, 1)
    SSSm <- measuresRSS(sample, i, 1000, stRSS) 
    USSm <- measures(sample, st.measure, i, 1000)
    SSSs <- measuresRSS(sample, i, 1000, markovRSS) 
    USSs <- measures(sample, markovInv, i, 1000)
    UnSSS <- measuresUnBiased(sample, i, 1000) 
    resultArray <- c(resultArray, c(SSSm[1]/1000, USSm[1]/1000, 
                                    SSSs[1]/1000, USSs[1]/1000, 
                                    UnSSS[1]/1000))
  }
  matrix(resultArray, 5, seqlen/step - 1)
}

step <- 5
simLoop(ex1, p, 5, 100) -> data

njFels <- read.csv("/Users/ataylor/research/colaboration/Tasmania/representation-theory-phylogenetics/NJFels.csv")
njbalanced <- read.csv("/Users/ataylor/research/colaboration/Tasmania/representation-theory-phylogenetics/NJBalanced.csv")
njFarris <- read.csv("/Users/ataylor/research/colaboration/Tasmania/representation-theory-phylogenetics/NJFarris.csv")
njFarrisStar <- read.csv("/Users/ataylor/research/colaboration/Tasmania/representation-theory-phylogenetics/NJStarFarris.csv")
njBalancedBig <- read.csv("/Users/ataylor/research/colaboration/Tasmania/representation-theory-phylogenetics/NJBalanceBig.csv")
njBalancedMid <- read.csv("/Users/ataylor/research/colaboration/Tasmania/representation-theory-phylogenetics/NJBalanceMid.csv")
njBalancedMid2 <- read.csv("/Users/ataylor/research/colaboration/Tasmania/representation-theory-phylogenetics/NJBalanceMid2.csv")

### Graphing the results.  Note that one needs to change only the first entry in 
### simLoop function above to use the different sets of parameters.  


#title <- "Variable Sequence Length - Star Tree"
plot(seq(10,100,step),data[1,], #main=title, 
     xlab="sequence length", ylab="percent correct",
     pch ="*", xlim=c(0,100), ylim=c(0,1), col="red")
par(new=T)
plot(seq(10,100,step),data[2,], #main=title, 
     xlab="sequence length", ylab="percent correct",  
     pch =".", xlim=c(0,100), ylim=c(0,1), col="red")
par(new=T)
plot(seq(10,100,step),data[3,], #main=title, 
     xlab="sequence length", ylab="percent correct",  
     pch ="*", xlim=c(0,100), ylim=c(0,1), col="blue")
par(new=T)
plot(seq(10,100,step),data[4,], #main=title, 
     xlab="sequence length", ylab="percent correct",  
     pch =".", xlim=c(0,100), ylim=c(0,1), col="blue")
par(new=T)
plot(seq(10,100,step),data[5,], #main=title, 
     xlab="sequence length", ylab="percent correct",  
     pch =".", xlim=c(0,100), ylim=c(0,1), col="purple")
par(new=T)
plot(njBalancedMid2$sequence.length,njBalancedMid2$NJ.Corrected/1000, #main=title, 
     xlab="sequence length", ylab="percent correct",  
     pch ="*", xlim=c(0,100), ylim=c(0,1), col="green")
lines(seq(10,100,step), data[1,], col="red", lty=1) #SSSminors
lines(seq(10,100,step), data[2,], col="red", lty=2) #USSminors
lines(seq(10,100,step), data[3,], col="blue", lty=1) #SSSsquangle
lines(seq(10,100,step), data[4,], col="blue", lty=2) #USSsquangle
lines(seq(10,100,step), data[5,], col="purple",lty=1) #Unbiased Signed Squangle
lines(njBalancedMid2$sequence.length,njBalancedMid2$NJ.Corrected/1000, col="green", lty = 1)#NJ
abline(a = 0.3333, b = 0, col = "orange", lty = 2)
#text(c(60,60,60),c(175, 100, 25),
#     labels = c("Green: Unbiased Signed Squangles",
#                "Red: Dashed Signed Minors/Solid Unsigned Minors",
#                "Blue: Dashed Signed Squangle/Solid Unsigned Squangle"))

par(new=F)

