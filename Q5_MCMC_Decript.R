require(R.matlab)

# text to decode
cvtext = c(19,15,16,24,18,17,16,18,16,8,1,19,20,26,15,16,6,9,14,2,16,2,18,21,16,19,3,16,18,13,1,19,14,16,18,3,
           2,16,15,26,27,16,6,14,9,6,10,17,16,24,27,1,27,16,17,15,1,19,10,19,3,20,16,15,26,19,1,15,27,27,3)
## Import Data ----
transition.matrix <- readMat('English_trans.mat')
# A = bigram transition A(i,j) = Prob(x_t = i| x_{t-1} = j)
A <- transition.matrix$A
# S = trigram transition A(i,j,k) = Prob(x_t = i| x_{t-1} = j, x_{t-2} = k)
S <- transition.matrix$S

## Functions -----
calc.log.plaus.bi <- function(A, f, cvtext){
  log.plaus <- 0
  for (k in 2:length(cvtext)){
    log.plaus = log.plaus + log(A[f[cvtext[k]],f[cvtext[k-1]]])
  }
  return(log.plaus)
}
calc.log.plaus.tri <- function(A, f, cvtext){
  log.plaus <- 0
  for (k in 3:length(cvtext)){
    log.plaus = log.plaus + log(S[f[cvtext[k]],f[cvtext[k-1]],f[cvtext[k-2]]])
  }
  return(log.plaus)
}
## Clean Data: remove zeros from the transition matrix and renormalize ----
tol <- 1e-6
# Clean A
colSums(A)
min(A)
A <- A + tol
for (j in 1:ncol(A)){
    A[,j] = A[,j]/colSums(A)[j]    
}
colSums(A)
min(A)

# Clean S
S = S + tol
for (k in 1:27) {
  for (l in 1:27) {
    sumprob <- 0
    for (m in 1:27) {
      sumprob <- sumprob + S[m,l,k];
    }
    if (sumprob>0) {
      for (m in 1:27) {
        S[m,l,k] <- S[m,l,k] / sumprob;    
      }
    }
  }  
}
## Settings ----
log.plaus.max <- -99999999
n <- 500000 # Number of iterations per chain
m <- 10 # Number of Chains to simulate
freq.report <- 20000 # Reporting frequency
plaus.hist <- matrix(data = 0, nrow = n / freq.report, ncol = m)
for (chains in 1:m) {
  cat('\n Chain # ', chains, '\n')
  # Initial function, picked at random
  f <- sample(x = c(1:27), size = 27, replace = FALSE)
  # Numbrer of iterations
  
  ## Loop ----
  for (counter in 1:n) {
    # Calculate Log likelihood of plausibility of function f
    log.plaus <- calc.log.plaus.bi(A, f, cvtext) #+ calc.log.plaus.bi(A, f, cvtext) 
    
    if (log.plaus > log.plaus.max) {
      log.plaus.max <- log.plaus
      f.max.plaus <- f
    }
    ## Make random transposition of  2 letters of f
    # Draw 2 random numbers
    change <- sample(1:27, 2)
    # Transpose those positions
    temp <- f[change[2]]
    f[change[2]] <- f[change[1]]
    f[change[1]] <- temp
    
    # Calculate new log likelihood of plausibility of new function f*
    log.plaus.2 <- calc.log.plaus.bi(A, f, cvtext) #+ calc.log.plaus.bi(A, f, cvtext)
    
    # If the plaus is worse, we flip a biased coin to see if we accept the new function. 
    if (log.plaus.2 < log.plaus) {
      coin <- runif(n = 1, min = 0, max = 1)
      if (coin > exp(x = log.plaus.2 ) / exp(x = log.plaus )) {
        # We switch back because we got tails
        temp <- f[change[2]]
        f[change[2]] <- f[change[1]]
        f[change[1]] <- temp
      }
    }
    if (counter %% freq.report == 0) {
      cat('Progress: ', (counter / n), ', Log Plausibility: ', log.plaus, 'Max: ',log.plaus.max, '\n')
      pos <- floor(counter / freq.report)
      plaus.hist[pos,chains] <- log.plaus
    }
    
  }
}
## Convert to text ----
ctext <- matrix(0, nrow = 1, ncol = length(cvtext))
# Apply the function f
for (k in 1:length(cvtext)) {
  ctext[k] <- f.max.plaus[cvtext[k]]
}
# Do the transformation
ctextconv <- matrix('0', nrow = 1, ncol = length(ctext))
for(i in 1:length(ctext)){
  if (ctext[i] == 27) {
    ctext[i] <- -64
  }  
  ctext[i] <- ctext[i] + 96
}
## Final output ----
cat('Best Log Plausibiilty: ', log.plaus.max, '\n')
ctextorig <- rawToChar(as.raw(cvtext+96))
ctextconv <- rawToChar(as.raw(ctext))
for (k in 1:length(ctextorig)) {
  cat(ctextorig[k])
}
cat('\n')
for (k in 1:length(ctextconv)) {
  cat(ctextconv[k])
}

## Plot ----
# get the range for the x and y axis 
yrange <- range(plaus.hist[1:(n/freq.report),]) 
# set up the plot 
plot(c(1, n/freq.report), yrange, type="n", xlab="Iteration (x20000)",
     ylab="Log Plausibility" ) 
linetype <- c(1:m) 
plotchar <- seq(18,18+m,1)
# add lines 
for (i in 1:m) { 
  lines(c(1:(n/freq.report)), plaus.hist[1:(n/freq.report),i], type="b", lwd=1.5,
        lty=linetype[i], pch=plotchar[i]) 
} 
# add a title and subtitle 
title("Results of Iterations, 10 different starting points, Log Plausibility Function")