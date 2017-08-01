# Here we continue studying the dependence of the convergence time on
# the variance, but now we consider s=1/2, 1 and 2, and functions
# cosx, cos(2x), cos(x/2), ..., cos(2^n x), cos(2^(-n) x)

# The number of functions is 2 * fn + 1
fn <- 3
fn21 <- fn * 2 + 1

# For function F(x) = 2^k * cos(2^k x) we have IF(x) = sin(2^k x), and
# PIF(s) = 2^k * exp(-2^(2k-1) * s^2)
ps2 = rep(0, fn21)
i <- 1
while (i <= fn21) {
  ps2[i] <- 2^(i - fn - 1)
  i <- i + 1
}

# The values of s we want to consider
Ss <- c( 0.3, 0.6, 1.2 )

# Epsilon
E <- 0.1

# Statistics we gather
# - the number of experiments for each s
N <- 10
# - one sample for each s
expXs <- list()
# - the resulting times
expTs <- list()
# - the resulting numbers of switching points
expSPNs <- list()

# Loop S
iS <- 1
while (iS <= length(Ss)) {

# the current standard deviation and the target values of integrals
s <- Ss[iS]
pifs <- c(fn)
i <- 1
while (i <= fn21) {
  k <- i - fn - 1
  pifs[i] <- ps2[i] * exp(- ps2[i]^2 * s^2 / 2)
  i <- i + 1
}

# initialize arrays to store experiment values for s there
Xs <- c()
Ts <- c()
expTs[[iS]] <- rep(0, N)
expSPNs[[iS]] <- rep(0, N)

# Loop N
iN <- 1
while (iN <= N) {

# Variables
iSP <- 1 # the index of the next switching point
theta <- 1 # the current theta
x <- 0 # the current x
t <- 0 # the current t
is <- rep(0, fn21) # the current values of integrals for each f

# Loop Zig-Zag
runZZ <- TRUE
while (runZZ || (iN == 1 && length(Ts) < 1000)) {

  # calculate the new value of x
  tau <- s * sqrt(-2 * log(runif(1)))
  xnew <- theta * tau

  # update the integrals
  i <- 1
  while (i <= fn21) {
    k <- i - fn - 1
    p2 <- ps2[i]
    is[i] <- is[i] + theta * (sin(p2 * xnew) - sin(p2 * x))
    i <- i + 1
  }

  # update t, including the time to move back to 0, and x
  t <- t + abs(x) + tau
  x <- xnew

  # if iN = 1 we collect the actual process
  if (iN == 1) {
    Xs <- c(Xs, x)
    Ts <- c(Ts, t)
  }

  # if the maximum measure is small enough we stop
  if (runZZ && max(abs(is / t / pifs - 1)) < E) {
    expTs[[iS]][iN] <- t
    expSPNs[[iS]][iN] <- iSP
    runZZ <- FALSE
  }

# The end of loop Zig-Zag
  iSP <- iSP + 1
  theta <- -theta
}

# The end of loop N
  iN <- iN + 1
}

# Obtain a sample from the process stored in one experiment
j <- 1
theta <- 1
n <- length(Ts)
expXs[[iS]] <- rep(0, n)
i <- 1
while (i <= n) {
  t <- i / (n + 1) * Ts[n]
  while (t > Ts[j]) { j <- j + 1; theta = -theta; }
  expXs[[iS]][i] <- Xs[j] - theta * (Ts[j] - t)
  i <- i + 1
}

# The end of loop S
  iS <- iS + 1
}

# Graphs
# Checking that we indeed generate N(0, s^2)
iS <- 1
while (iS <= length(Ss)) {
  plot(density(expXs[[iS]], bw = "SJ"))
  normalxs <- seq(min(expXs[[iS]]) - 1, max(expXs[[iS]]) + 1, length = 100)
  normalys <- exp(-1 / 2 / Ss[iS]^2 * normalxs^2) / sqrt(2 * pi) / Ss[iS]
  lines(normalxs, normalys, col = "red", lwd = 2)
  iS<- iS + 1
}
# Plotting average T against s
plot(Ss, sapply(expTs, mean), type = "l", col = "blue")

# Report results
N
rbind(Ss, sapply(expSPNs, mean), sapply(expSPNs, sd))
rbind(Ss, sapply(expTs, mean), sapply(expTs, sd))
