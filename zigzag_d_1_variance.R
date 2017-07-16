# Finctions to use to evaluate convergence
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

# The values of s we want to consider
Ss <- c( 0.125, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5 )

# Epsilon
E <- 0.05

# Statistics we gather
# - the number of experiments for each s
N <- 1000
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
pifs <- c(length(Fs))
i <- 1
while (i <= length(Fs)) { pifs[i] <- PIFs[[i]](s); i <- i + 1; }

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
is <- rep(0, length(Fs)) # the current values of integrals for each f

# Loop Zig-Zag
runZZ <- TRUE
while (runZZ || (iN == 1 && length(Ts) < 1000)) {

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
