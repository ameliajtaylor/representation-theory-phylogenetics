### This file is written assuming we generate data for the 12.34 tree.
### Implements both residual sum of squares and not. 
### This require the function unBiasedSQ so we need to source the following file.  This 
### path is set to match where the file is found in my directory structure. 

source("/Users/ataylor/research/colaboration/Tasmania/representation-theory-phylogenetics/unbiasedSquangleSquared-Jez.R")
options(expressions = 10000)

# Function for producing 2x2 markov matrices.
mark <- function(a,b) matrix(c(1-a,a,b,1-b),2,2)

# Function for producing two-taxa tree.

taxa2 <- function(a1,b1,a2,b2,p) {
  
  D <- matrix(c(p,0,0,1-p),2,2)
  M1 <- mark(a1,b1)
  M2 <- mark(a2,b2)
  ans <- M1%*%D%*%t(M2)
  
  ans
  
}

# Function for branching into 12/34 and flattening at the same time.

branch12.34 <- function(P) {
  
  ans <- matrix(0,4,4)
  ans[c(1,4),c(1,4)] <- P
  ans
  
}

# Final function. Inital root distribution is (p,(1-p))

tree12.34 <- function(a5,b5,a1,b1,a2,b2,a3,b3,a4,b4,p){
  
  # Root the tree at 12 divergence.
  P <- taxa2(a5,b5,0,0,p)
  P <- branch12.34(P)
  # Leaf action
  M12 <- mark(a1,b1)%x%mark(a2,b2)
  M34 <- mark(a3,b3)%x%mark(a4,b4)
  
  P <- M12%*%P%*%t(M34)
  
  P
  
}

# Multinomial sampling for sequence length N:

F12.34 <- function(par,p,N,R,t = 1){
  
  P <- tree12.34(par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8],par[9],par[10],p)
  if (t == 2) P <- convert(P, c(1, 9, 3, 11, 5, 13, 7, 15, 2, 10, 4, 12, 6, 14, 8, 16))
  if (t == 3) P <- convert(P,c(1, 5, 3, 7, 2, 6, 4, 8, 9, 13, 11, 15, 10, 14, 12, 16))
  F <- rmultinom(R,N,P)  #rmultinom treats P as a vector
  #  F <- matrix(F,4,4)
  F
}

######  Above we are building trees and samples
######  Below we are constructing measures

# Takes values from the two "incorrect" flattenings and produces the 
#residual sum of squares. This is the decision rule. 
rss <- function(x,y){
  if ( (x -y) > 0) {1/2*(x + y)^2} 
  else (x^2+y^2)
}



# Computes Markov Invariants by converting to the new
# basis via hadamard and computing the appropriate minor.

markovInv <- function(P){
  cB <- solve(matrix(c(1/2,-1/2,1/2,1/2),2,2))
  newBasis <- (cB%x%cB)%*%P%*%t(cB%x%cB)
  det(newBasis[-1,-1])
}
# (14)P, I think --- c(1, 9, 3, 11, 5, 13, 7, 15, 2, 10, 4, 12, 6, 14, 8, 16)
# (23)P --- c(1, 2, 5, 6, 3, 4, 7, 8, 9, 10, 13, 14, 11, 12, 15, 16)
# Takes a matrix P, assumes it is T1 and gives the markov invariants Q1, Q2, Q3.  
# Sanity Check:  The sum of the vector entries should be zero in theory. 
markovInv.vec <- function(P) {
  c(markovInv(P),
    -1*markovInv(convert(P,c(1, 9, 3, 11, 5, 13, 7, 15, 2, 10, 4, 12, 6, 14, 8, 16))), 
    -1*markovInv(convert(P,c(1, 5, 3, 7, 2, 6, 4, 8, 9, 13, 11, 15, 10, 14, 12, 16))) )
}

# Takes P, produces the Q1, Q2, Q3 from markInv.vec and then does the decision process. 
# returns a vector with 3 elements.  The tree with the smallest residual, that residual, and the
# difference needed for our test statistic, allowing both tree counting and the use 
# of the statistc.
  
markovRSS <- function(P) {
  Q <- markovInv.vec(P)
  rs <- c(rss(Q[3], Q[2]), rss(Q[1], Q[3]), rss(Q[2], Q[1]))
  #Comment out these lower 3 and uncomment last one for Variance runs
  m = which.min(rs)
  if (length(which(rs == rs[m])) > 1) {c(4, rs[m], rs[m]-rs[1])}
  else c(m, rs[m], min(rs[2],rs[3]) - rs[1],rs)
  #rs
}

markovUnBiased <- function(P) {
  T2 <- convert(P,c(1, 9, 3, 11, 5, 13, 7, 15, 2, 10, 4, 12, 6, 14, 8, 16))
  T3 <- convert(P,c(1, 5, 3, 7, 2, 6, 4, 8, 9, 13, 11, 15, 10, 14, 12, 16))
  Q <- markovInv.vec(P)
  u <- Q[3] - Q[2]
  v <- Q[1] - Q[3]
  w <- Q[2] - Q[1]
  s1 <- unBiasedSq(P)
  s2 <- unBiasedSq(T2)
  s3 <- unBiasedSq(T3)
  rs <-c(0,0,0)
  if (u > 0) {
    if (v < 0) {
      if (w < 0) {rs[1] <- 1}
      else{ if(s1 < s3) {rs[1] <- 1}
            else{rs[3] <- 1}
      }
    }
    if(s1 < s2){rs[1] <- 1}
    else{rs[2] <- 1}
    }
  else{
    if (v > 0){
      if (w < 0){rs[2] <- 1}
      else {
        if(s2 < s3){rs[2] <- 1}
        else{rs[3] <- 1}
      }
    }
    else{rs[3]<-1}
  }
  rs
}


# generate the usual minors from a matrix P.  
# Input is a matrix of probabilites. 

minors <- function(P){
  m <- matrix(c(1:16), 4,4)
  for(i in 1:4){
    for(j in 1:4){
      m[i,j]<-det(P[-i,-j])
    }
  }
  m
}

# Computes the residual sum of squares for the minors from the matrix.  

stRSS <- function(P) {
  
  # S1 and S2 set the minors to be an expected sign so that it is clear what order 
  # to place the minors into the RSS decision function rss.  
  S1 <- matrix(c(1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1), 4, 4)
  S2 <- matrix(c(-1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1), 4, 4)
  m1 <- minors(P)
  m2 <- minors(convert(P,c(1, 9, 3, 11, 5, 13, 7, 15, 2, 10, 4, 12, 6, 14, 8, 16))) 
  m3 <- minors(convert(P,c(1, 5, 3, 7, 2, 6, 4, 8, 9, 13, 11, 15, 10, 14, 12, 16)))
  
  # Set up so that q2 is all positive and q3 is all negative. 
  q2 <- S1*m2
  q3 <- S1*m3
  Q1 <- sum(c(rss(q2[1,1], q3[1,1]),
              rss(q2[1,3], q3[1,2]),
              rss(q2[3,1], q3[3,1]),
              rss(q2[3,3], q3[3,2]),
              rss(q2[1,2], q3[2,1]),
              rss(q2[1,4], q3[2,2]),
              rss(q2[3,2], q3[4,1]),
              rss(q2[3,4], q3[4,2]),
              rss(q2[2,1], q3[1,3]),
              rss(q2[2,3], q3[1,4]), 
              rss(q2[4,1], q3[3,3]),
              rss(q2[4,3], q3[3,4]),
              rss(q2[2,2], q3[2,3]),
              rss(q2[2,4], q3[2,4]),
              rss(q2[4,2], q3[4,3]),
              rss(q2[4,4], q3[4,4])))
  
  # Set up so that q1 is all positive and q2 is all negative.
  q1 <- S1*m1
  q3 <- S2*m3
  Q2 <- sum(c(rss(q1[1,1], q3[1,1]),
              rss(q1[2,1], q3[1,2]),
              rss(q1[3,1], q3[3,1]),
              rss(q1[4,1], q3[3,2]),
              rss(q1[1,2], q3[2,1]),
              rss(q1[2,2], q3[2,2]),
              rss(q1[3,2], q3[4,1]),
              rss(q1[4,2], q3[4,2]),
              rss(q1[1,3], q3[1,3]),
              rss(q1[2,3], q3[1,4]),
              rss(q1[3,3], q3[3,3]),
              rss(q1[4,3], q3[3,4]),
              rss(q1[1,4], q3[2,3]),
              rss(q1[2,4], q3[2,4]),
              rss(q1[3,4], q3[4,3]),
              rss(q1[4,4], q3[4,4])))
  
  # Set up so that q1 is all positive and q2 is all negative.
  q1 <- S2*m1
  q2 <- S1*m2
  Q3 <- sum(c(rss(q1[1,1], q2[1,1]),
              rss(q1[2,1], q2[1,3]),
              rss(q1[3,1], q2[3,1]),
              rss(q1[4,1], q2[3,3]),
              rss(q1[1,2], q2[1,2]),
              rss(q1[2,2], q2[1,4]),
              rss(q1[3,2], q2[3,2]),
              rss(q1[4,2], q2[3,4]),
              rss(q1[1,3], q2[2,1]),
              rss(q1[2,3], q2[2,3]),
              rss(q1[3,3], q2[4,1]),
              rss(q1[4,3], q2[4,3]),
              rss(q1[1,4], q2[2,2]),
              rss(q1[2,4], q2[2,4]),
              rss(q1[3,4], q2[4,2]),
              rss(q1[4,4], q2[4,4])))
  rs <- c(Q1, Q2, Q3)
  
  # Next two commands are set up to check we don't get bias from the way R 
  # computes minimums and records ties as a fourth unidentifiable tree. 
  ### Comment out below and uncomment last one for running VarianceRSS. 
  m <- which.min(rs)
  if (length(which(rs == rs[m])) > 1) {c(4, rs[m], rs[m]-rs[1])}
  
  # Returns the number of the tree with smallest RSS, the RSS for that tree 
  # and the computation for the numerator of the test statistic. 
  
  else c(m, rs[m], min(rs[2],rs[3]) - rs[1], rs)
  #rs
}

### Same as above, but only uses the minors we hypothesize to have biggest
### variance. These generate the 2 dimensional subspace under the leaf
### permutations.

bigVarRSS <- function(P) {
  
  # S1 and S2 set the minors to be an expected sign so that it is clear what order 
  # to place the minors into the RSS decision function rss.  
  S1 <- matrix(c(1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1), 4, 4)
  S2 <- matrix(c(-1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1), 4, 4)
  m1 <- minors(P)
  m2 <- minors(convert(P,c(1, 9, 3, 11, 5, 13, 7, 15, 2, 10, 4, 12, 6, 14, 8, 16))) 
  m3 <- minors(convert(P,c(1, 5, 3, 7, 2, 6, 4, 8, 9, 13, 11, 15, 10, 14, 12, 16)))
  
  # Set up so that q2 is all positive and q3 is all negative. 
  q2 <- S1*m2
  q3 <- S1*m3
  Q1 <- sum(c(rss(q2[3,3], q3[3,2]),
              rss(q2[2,2], q3[2,3])))
  
  # Set up so that q1 is all positive and q2 is all negative.
  q1 <- S1*m1
  q3 <- S2*m3
  Q2 <- sum(c(rss(q1[2,2], q3[2,2]),
              rss(q1[3,3], q3[3,3])))
  
  # Set up so that q1 is all positive and q2 is all negative.
  q1 <- S2*m1
  q2 <- S1*m2
  Q3 <- sum(c(rss(q1[3,2], q2[3,2]),
              rss(q1[2,3], q2[2,3])))
  rs <- c(Q1, Q2, Q3)
  
  # Next two commands are set up to check we don't get bias from the way R 
  # computes minimums and records ties as a fourth unidentifiable tree. 
  m <- which.min(rs)
  if (length(which(rs == rs[m])) > 1) {c(4, rs[m], rs[m]-rs[1])}
  
  # Returns the number of the tree with smallest RSS, the RSS for that tree 
  # and the computation for the numerator of the test statistic. 
  else c(m, rs[m], min(rs[2],rs[3]) - rs[1])
}


### Same as above, but only uses the minors we hypothesize to have the
### smallest variance. These generate one of the four dimensional
### subspaces under the leaf permutations.

smallVarRSS <- function(P) {
  
  # S1 and S2 set the minors to be an expected sign so that it is clear what order 
  # to place the minors into the RSS decision function rss.  
  S1 <- matrix(c(1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1), 4, 4)
  S2 <- matrix(c(-1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1), 4, 4)
  m1 <- minors(P)
  m2 <- minors(convert(P,c(1, 9, 3, 11, 5, 13, 7, 15, 2, 10, 4, 12, 6, 14, 8, 16))) 
  m3 <- minors(convert(P,c(1, 5, 3, 7, 2, 6, 4, 8, 9, 13, 11, 15, 10, 14, 12, 16)))
  
  # Set up so that q2 is all positive and q3 is all negative. 
  q2 <- S1*m2
  q3 <- S1*m3
  Q1 <- sum(c(rss(q2[1,4], q3[2,2]),
              rss(q2[3,2], q3[4,1]),
              rss(q2[2,3], q3[1,4]), 
              rss(q2[4,1], q3[3,3])))
  
  # Set up so that q1 is all positive and q2 is all negative.
  q1 <- S1*m1
  q3 <- S2*m3
  Q2 <- sum(c(rss(q1[4,1], q3[3,2]),
              rss(q1[3,2], q3[4,1]),
              rss(q1[2,3], q3[1,4]),
              rss(q1[1,4], q3[2,3])))
  
  # Set up so that q1 is all positive and q2 is all negative.
  q1 <- S2*m1
  q2 <- S1*m2
  Q3 <- sum(c(rss(q1[4,1], q2[3,3]),
              rss(q1[2,2], q2[1,4]),
              rss(q1[3,3], q2[4,1]),
              rss(q1[1,4], q2[2,2])))
  
  rs <- c(Q1, Q2, Q3)
  
  # Next two commands are set up to check we don't get bias from the way R 
  # computes minimums and records ties as a fourth unidentifiable tree. 
  m <- which.min(rs)
  if (length(which(rs == rs[m])) > 1) {c(4, rs[m], rs[m]-rs[1])}

  # Returns the number of the tree with smallest RSS, the RSS for that tree 
  # and the computation for the numerator of the test statistic. 
  else c(m, rs[m], min(rs[2],rs[3]) - rs[1])
}

### Same as above, but only uses the minors from the subspaces we think do
### not contain the ID and SGN reps.

noIdSgnRSS <- function(P) {
  
  # S1 and S2 set the minors to be an expected sign so that it is clear what order 
  # to place the minors into the RSS decision function rss.  
  S1 <- matrix(c(1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1), 4, 4)
  S2 <- matrix(c(-1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1), 4, 4)
  m1 <- minors(P)
  m2 <- minors(convert(P,c(1, 9, 3, 11, 5, 13, 7, 15, 2, 10, 4, 12, 6, 14, 8, 16))) 
  m3 <- minors(convert(P,c(1, 5, 3, 7, 2, 6, 4, 8, 9, 13, 11, 15, 10, 14, 12, 16)))
  
  # Set up so that q2 is all positive and q3 is all negative. 
  q2 <- S1*m2
  q3 <- S1*m3
  Q1 <- sum(c(rss(q2[1,1], q3[1,1]),
              rss(q2[1,3], q3[1,2]),
              rss(q2[3,1], q3[3,1]),
              rss(q2[1,2], q3[2,1]),
              rss(q2[3,4], q3[4,2]),
              rss(q2[2,1], q3[1,3]),
              rss(q2[4,3], q3[3,4]),             
              rss(q2[2,4], q3[2,4]),
              rss(q2[4,2], q3[4,3]),
              rss(q2[4,4], q3[4,4])))
  
  # Set up so that q1 is all positive and q2 is all negative.
  q1 <- S1*m1
  q3 <- S2*m3
  Q2 <- sum(c(rss(q1[1,1], q3[1,1]),
              rss(q1[2,1], q3[1,2]),
              rss(q1[3,1], q3[3,1]),
              rss(q1[1,2], q3[2,1]),
              rss(q1[4,2], q3[4,2]),
              rss(q1[1,3], q3[1,3]),
              rss(q1[4,3], q3[3,4]),
              rss(q1[2,4], q3[2,4]),
              rss(q1[3,4], q3[4,3]),
              rss(q1[4,4], q3[4,4])))
  
  # Set up so that q1 is all positive and q2 is all negative.
  q1 <- S2*m1
  q2 <- S1*m2
  Q3 <- sum(c(rss(q1[1,1], q2[1,1]),
              rss(q1[2,1], q2[1,3]),
              rss(q1[3,1], q2[3,1]),
              rss(q1[1,2], q2[1,2]),
              rss(q1[4,2], q2[3,4]),
              rss(q1[1,3], q2[2,1]),
              rss(q1[4,3], q2[4,3]),
              rss(q1[2,4], q2[2,4]),
              rss(q1[3,4], q2[4,2]),
              rss(q1[4,4], q2[4,4])))
  rs <- c(Q1, Q2, Q3)
  
  # Next two commands are set up to check we don't get bias from the way R 
  # computes minimums and records ties as a fourth unidentifiable tree. 
  m <- which.min(rs)
  if (length(which(rs == rs[m])) > 1) {c(4, rs[m], rs[m]-rs[1])}
  
  # Returns the number of the tree with smallest RSS, the RSS for that tree 
  # and the computation for the numerator of the test statistic. 
  else c(m, rs[m], min(rs[2],rs[3]) - rs[1])
}


# A measure using the standard basis for minors. 
# input is P

st.measure <- function(P){
  M <- minors(P)
#  sqrt(sum(M^2))
  1/2*(sum(M^2))
}


#convert measures between trees 12.34, 13.24, and 14.23
# Accepts a matrix of probabilites or counts and a permutation vector.  
convert <- function(P,perm){
  dim(P) <- NULL
  Pnew <- P[perm]
  dim(Pnew) <- c(4,4)
  Pnew
}

# If f is st.measure then this gives the sum of squares for the
# standard minors.  Allows other measures given by f. We could
# then get the sum of the squares for the Markov invariant using
# the appropriate function f.

measures <- function(P,f,N,R) { 
  
  #initialize vectors of measures for the 3 flattenings and the corresponding testStat
  d1 <- NULL  
  d2 <- NULL
  d3 <- NULL
  tS <- NULL
  count <- c(0,0,0,0)
  
  #The loop through all the samples getting the residual sum of squares, 
  # counting and computing the key info for the test statistic. 
  for (i in 1:R){
    P12.34 <- matrix(P[,i],4,4)
    #convert P, known to be on 12.34 to 13.24 and 14.23 flattenings
    P13.24 <- convert(P12.34,c(1, 9, 3, 11, 5, 13, 7, 15, 2, 10, 4, 12, 6, 14, 8, 16))  
    P14.23 <- convert(P12.34,c(1, 5, 3, 7, 2, 6, 4, 8, 9, 13, 11, 15, 10, 14, 12, 16))
    d1[i] <- f(P12.34)
    d2[i] <- f(P13.24)
    d3[i] <- f(P14.23)
    tS[i] <- min(d2[i],d3[i])-d1[i]
    delta <- c(d1[i],d2[i],d3[i])
    m = which.min(abs(delta))
    if (length(which(delta == delta[m])) > 1) {count[4] <- count[4]+1}
    else count[m] <- count[m]+1
  }
  c(count, mean(tS)/sd(tS))           
}

# Using signed least squares --- particular function given by f --- counts trees and returns 
# the test statistic from a sample run using F12.34
# standard minors stRSS
# markov invariant markovRSS
# small/big variance smallVarRSS bigVarRSS
# the modules disjoint from those that contain Id and Sgn --- medium variance. 

measuresRSS <- function(P,N,R,f) { 
  #  d1 <- NULL  
  count <- c(0,0,0,0)
  tS <- NULL
  
  #The loop through all the samples getting the residual sum of squares, 
  # counting and computing the key info for the test statistic. 
  for (i in 1:R){
    P12.34 <- matrix(P[,i],4,4)/N
    #convert P, known to be on 12.34 to 13.24 and 14.23 flattenings
    m <- f(P12.34)
    count[m[1]] <- count[m[1]]+1
    tS[i] <- m[3] 
  }
  c(count, mean(tS)/sd(tS))             
} 

measuresUnBiased <- function(P,N,R) { 
  #  d1 <- NULL  
  count <- c(0,0,0,0)
  tS <- NULL
  
  #The loop through all the samples getting the residual sum of squares, 
  # counting and computing the key info for the test statistic. 
  for (i in 1:R){
    P12.34 <- matrix(P[,i],4,4)
    #convert P, known to be on 12.34 to 13.24 and 14.23 flattenings
    m <- markovUnBiased(P12.34)
    test <- which(m == 1)
    count[test] <- count[test]+1
  }        
  count
} 