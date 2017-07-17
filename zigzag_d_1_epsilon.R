# In this experiment we consider simulating N(0,s^2)
# until the measure of convergence becomes less than epsilon (E)
# We do this for different values of E

# The standard deviation of the target normal distribution
s <- 1

# Functions to use to evaluate convergence
Fs <- c(
  function (x) { abs(x); }
  , function (x) { x^2; }
  , function (x) { 1 / sqrt(abs(x)); }
#  , function (x) { cos(x); }
)

# Explicit integrals of these functions
# We will need int [t from 0 to tau] f(x0 + theta t) dt =
#              int [x from x0 to x0 + theta T] f(x) / theta dx =
#              theta * (F(xnew) - F(x0)),
# where F(x) = int f(x) dx is provided here for each f
# (in the future we may need to integrate numerically instead
# in more complicated cases)
IFs <- c(
  function (x) { x * abs(x) / 2; }
  , function (x) { x^3 / 3; }
  , function (x) { 2 * sign(x) * sqrt(abs(x)); }
  , function (x) { sin(x); }
)

# pi(f) for each f w.r.t. N(0, s^2)
PIFs <- c(
  function (s) { sqrt(2 / pi) * s; }
  , function (s) { s^2; }
  , function (s) { 1.720079974649039070752407248933115962110249350333550144551 / sqrt(s); }
  , function (s) { exp(-s^2 / 2); }
)

# Epsilon values
# These should be sorted from lowest to largest
Es <- c(seq(0.01, 0.1, length = 10), seq(0.15, 0.5, length = 8))

# Statistics we gather
# - the number of experiments for each E
N <- 1000
# - one sample
expXs <- c()
# - the resulting times
expTs <- list()
# - the resulting numbers of switching points
expSPNs <- list()

# the target values of integrals
pifs <- c(length(Fs))
i <- 1
while (i <= length(Fs)) { pifs[i] <- PIFs[[i]](s); i <- i + 1; }

# Initialize arrays to store experiment values
# these will store one process to see that the algorithm works
Xs <- c()
Ts <- c()
# these will store the resulting times and switching points for each E
iE <- 1
while (iE <= length(Es)) {
  expTs[[iE]] <- rep(0, N)
  expSPNs[[iE]] <- rep(0, N)
  iE <- iE + 1
}

# Loop N
iN <- 1
while (iN <= N) {

# Variables
iSP <- 1 # the index of the next switching point
theta <- 1 # the current theta
x <- 0 # the current x
t <- 0 # the current t
is <- rep(0, length(Fs)) # the current values of integrals for each f

# the current standard deviation
iE <- length(Es)
E <- Es[iE]

# Loop Zig-Zag
while (iE > 0 || (iN == 1 && length(Ts) < 1000)) {

  # calculate the new value of x
  tau <- s * sqrt(-2 * log(runif(1)))
  xnew <- theta * tau

  # update the integrals
  i <- 1
  while (i <= length(Fs)) {
    is[i] <- is[i] + theta * (IFs[[i]](xnew) - IFs[[i]](x))
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
  maxmeasure <- max(abs(is / t / pifs - 1))
  while (iE > 0 && maxmeasure < E) {
    expTs[[iE]][iN] <- t
    expSPNs[[iE]][iN] <- iSP
    iE <- iE - 1
    E <- Es[iE]
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
expXs <- rep(0, n)
i <- 1
while (i <= n) {
  t <- i / (n + 1) * Ts[n]
  while (t > Ts[j]) { j <- j + 1; theta = -theta; }
  expXs[i] <- Xs[j] - theta * (Ts[j] - t)
  i <- i + 1
}

# Graphs
# Checking that we indeed generate N(0, 1)
plot(density(expXs, bw = "SJ"))
normalxs <- seq(min(expXs) - 1, max(expXs) + 1, length = 100)
normalys <- exp(-1 / 2 / s^2 * normalxs^2) / sqrt(2 * pi) / s
lines(normalxs, normalys, col = "red", lwd = 2)
# Plotting average T against E
plot(Es, sapply(expTs, mean), type = "l", col = "blue")

# Report results
N
rbind(Es, sapply(expSPNs, mean), sapply(expSPNs, sd))
rbind(Es, sapply(expTs, mean), sapply(expTs, sd))
Es * sapply(expTs, mean)
