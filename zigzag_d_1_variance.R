# Finctions to use to evaluate convergence
Fs <- c(
  function (x) { abs(x); },
  function (x) { x^2; },
  function (x) { abs(x)^3; },
  function (x) { 1 / x^2; },
  function (x) { cos(x); }
)

# Explicit integrals of these functions
# We will need int [t from 0 to tau] f(x0 + theta t) dt =
#              int [x from x0 to x0 + theta T] f(x) / theta dx =
#              theta * (F(xnew) - F(x0)),
# where F(x) = int f(x) dx is provided here for each f
# (in the future we may need to integrate numerically instead
# in more complicated cases)
IFs <- c(
  function (x) { x * abs(x) / 2; },
  function (x) { x^3 / 3; },
  function (x) { x * abs(x)^3 / 4; },
  function (x) { -1 / x; },
  function (x) { sin(x); }
)

# pi(f) for each f w.r.t. N(0, s^2)
PIFs <- c(
  function (s) { sqrt(2 / pi) * s; },
  function (s) { s^2; },
  function (s) { sqrt(8 / pi) * s^3; },
  function (s) { -1 / s^2; },
  function (s) { exp(-s^2 / 2); }
)

# The values of s we want to consider
Ss <- c( 0.5, 1, 2 )

# Epsilon
E <- 0.01

# Statistics we gather
# - the number of experiments for each s
N <- 1000
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

# initialize expTs[[iS]] and expSPNs[[iS]] to store values there
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
while (TRUE) {

  # generate tau
  tau <- s * sqrt(-2 * log(runif(1)))

  # calculate the new value of x
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

  # if the maximum measure is small enough we stop
  if (max(is / t / pifs - 1) < E) {
    expTs[[iS]][iN] <- t
    expSPNs[[iS]][iN] <- iSP
    break
  }

# The end of loop Zig-Zag
  iSP <- iSP + 1
  theta <- -theta
}

# The end of loop N
  iN <- iN + 1
}

# The end of loop S
  iS <- iS + 1
}

# Densities
#iS <- 1
#while (iS <= length(Ss)) {
#  plot(density(expTs[[iS]], bw = "SJ"))
#  iS <- iS + 1
#}

# Means and standard deviations
rbind(Ss, sapply(expTs, mean), sapply(expTs, sd))
