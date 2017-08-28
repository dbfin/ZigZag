# Here we study the convergence of the Zig-Zag algorithm for N((0,0), (1, a, a, s^2))

# The values of s to consider
Ss <- c(1, 1.5)

# The values of a to consider
As <- c(0, 0.25, 0.5)

# The number of experiments for each s and a
N <- 100

# Finctions f to use to evaluate convergence
Fs <- c(
  function (x, y) { abs(x); },
  function (x, y) { x^2; },
  function (x, y) { 1 / sqrt(abs(x)); },
  function (x, y) { abs(y); },
  function (x, y) { y^2; },
  function (x, y) { 1 / sqrt(abs(y)); },
  function (x, y) { abs(x * y); }
)

# Explicit integrals of these functions
# We will need int [t from 0 to tau] f(x0 + theta_x t, y0 + theta_y t) dt
# When f depends on, say, x only, this becomes
#              int [t from 0 to tau] f(x0 + theta t) dt =
#              int [x from x0 to x1] f(x) / theta dx =
#              theta * (F(x1) - F(x0)),
# where F(x) = int f(x) dx
IFs <- c(
  function (x0, y0, x1, y1, theta_x, theta_y, tau) { theta_x * (x1 * abs(x1) - x0 * abs(x0)) / 2; },
  function (x0, y0, x1, y1, theta_x, theta_y, tau) { theta_x * (x1^3 - x0^3) / 3; },
  function (x0, y0, x1, y1, theta_x, theta_y, tau) { theta_x * (sign(x1) * sqrt(abs(x1)) - sign(x0) * sqrt(abs(x0))) * 2; },
  function (x0, y0, x1, y1, theta_x, theta_y, tau) { theta_y * (y1 * abs(y1) - y0 * abs(y0)) / 2; },
  function (x0, y0, x1, y1, theta_x, theta_y, tau) { theta_y * (y1^3 - y0^3) / 3; },
  function (x0, y0, x1, y1, theta_x, theta_y, tau) { theta_y * (sign(y1) * sqrt(abs(y1)) - sign(y0) * sqrt(abs(y0))) * 2; },
  function (x0, y0, x1, y1, theta_x, theta_y, tau) {
    boundaries <- sort(c(-x0 / theta_x, -y0 / theta_y))
    boundaries <- boundaries[which(boundaries > 0 & boundaries < tau)]
    boundaries <- c(0, boundaries, tau)
    i <- 1
    t0 <- 0
    if (x0 == 0 || y0 == 0) {
      s <- sign((x0 + theta_x * boundaries[2]) * (y0 + theta_y * boundaries[2]))
    }
    else {
      s <- sign(x0 * y0)
    }
    integral <- 0
    while (i < length(boundaries)) {
      t1 <- boundaries[i + 1]
      integral <- integral + s * (x0 * y0 * (t1 - t0) + (theta_x * y0 + theta_y * x0) * (t1^2 - t0^2) / 2 + theta_x * theta_y * (t1^3 - t0^3) / 3)
      s <- -s
      t0 <- t1
      i <- i + 1
    }
    integral
  }
)

# pi(f)
PIFs <- c(
  function (s, a) { sqrt(2 / pi); },
  function (s, a) { 1; },
  function (s, a) { 1.720079974649039070752407248933115962110249350333550144551; },
  function (s, a) { sqrt(2 / pi) * s; },
  function (s, a) { s^2; },
  function (s, a) { 1.720079974649039070752407248933115962110249350333550144551 / sqrt(s); },
  function (s, a) { 2 * s / pi * (a * asin(a) + sqrt(1 - a^2)); }
)

# Epsilon
E <- 0.05

# Statistics we gather
# - one sample for each s and a
expXs <- array(dim = c(length(Ss), length(As), N))
expYs <- array(dim = c(length(Ss), length(As), N))
# - the resulting times
expTs <- array(dim = c(length(Ss), length(As), N))
# - the resulting numbers of switching points
expSPNs <- array(dim = c(length(Ss), length(As), N))

proc_time <- proc.time()

# Loop S
iS <- 1
while (iS <= length(Ss)) {

# Loop A
iA <- 1
while (iA <= length(As)) {

  # the current s and a, and the target values of integrals
  s <- Ss[iS]
  a <- As[iA]
  pifs <- c(length(Fs))
  i <- 1
  while (i <= length(Fs)) { pifs[i] <- PIFs[[i]](s, a); i <- i + 1; }

  # initialize arrays to store experiment values for s and a there
  Xs <- rep(0, 4 * N)
  Ys <- rep(0, 4 * N)
  Ts <- rep(0, 4 * N)

  # Loop N
  iN <- 1
  while (iN <= N) {
    
    # Variables
    iSP <- 1 # the index of the next switching point
    theta_x <- 1 # the current theta_x
    theta_y <- 1 # the current theta_y
    x <- 0 # the current x
    y <- 0 # the current y
    t <- 0 # the current t
    is <- rep(0, length(Fs)) # the current values of integrals for each f

    # Constants
    ca <- a / s
    cb <- 1 / s^2
    cc <- 1 / (1 - a^2)
    
    # Loop Zig-Zag
    runZZ <- TRUE
    while (runZZ || (iN == 1 && iSP <= 4 * N)) {

      # calculate the new values of x and y
      p_x <- cc * theta_x * (x - ca * y)
      p_y <- cc * theta_y * (cb * y - ca * x)
      q_x <- cc * (1 - ca * theta_x * theta_y)
      q_y <- cc * (cb - ca * theta_x * theta_y)
      tau_x <- -1
      if (q_x > 0) {
        tau_x <- max(0, -p_x / q_x)
        p_x <- max(0, p_x)
        tau_x <- tau_x + (sqrt(p_x^2 - 2 * q_x * log(runif(1))) - p_x) / q_x
      }
      else if (q_x < 0) {
        if (p_x > 0) {
          u <- runif(1)
          if (u >= exp(p_x^2 / q_x / 2)) {
            tau_x <- p_x^2 - 2 * q_x * log(u)
          }
        }
      }
      else {
        if (p_x > 0) {
          tau_x <- -log(runif(1)) / p_x
        }
      }
      tau_y <- -1
      if (q_y > 0) {
        tau_y <- max(0, -p_y / q_y)
        p_y <- max(0, p_y)
        tau_y <- tau_y + (sqrt(p_y^2 - 2 * q_y * log(runif(1))) - p_y) / q_y
      }
      else if (q_y < 0) {
        if (p_y > 0) {
          u <- runif(1)
          if (u >= exp(p_y^2 / q_y / 2)) {
            tau_y <- p_y^2 - 2 * q_y * log(u)
          }
        }
      }
      else {
        if (p_y > 0) {
          tau_y <- -log(runif(1)) / p_y
        }
      }
      if (tau_x < 0) {
        tau <- tau_y
      }
      else if (tau_y < 0) {
        tau <- tau_x
      }
      else {
        tau <- min(tau_x, tau_y)
      }
      x1 <- x + theta_x * tau
      y1 <- y + theta_y * tau

      # update the integrals
      i <- 1
      while (i <= length(Fs)) {
        is[i] <- is[i] + IFs[[i]](x, y, x1, y1, theta_x, theta_y, tau)
        i <- i + 1
      }
      
      # update t, x and y
      t <- t + tau
      x <- x1
      y <- y1
      
      # if iN = 1 we collect the actual process
      if (iN == 1 && iSP <= 4 * N) {
        Xs[iSP] <- x
        Ys[iSP] <- y
        Ts[iSP] <- t
      }
      
      # if the maximum measure is small enough we stop
      if (runZZ && max(abs(is / t / pifs - 1)) < E) {
        expTs[iS, iA, iN] <- t
        expSPNs[iS, iA, iN] <- iSP
        runZZ <- FALSE
      }
      
      # The end of loop Zig-Zag
      iSP <- iSP + 1
      if (tau == tau_x) { theta_x <- -theta_x; }
      else { theta_y <- -theta_y; }
    }
    
    # The end of loop N
    iN <- iN + 1
  }
  
  # Obtain a sample from the process stored in one experiment
  tj_1 <- 0
  xj_1 <- 0
  yj_1 <- 0
  j <- 1
  n <- length(Ts)
  i <- 1
  while (i <= N) {
    t <- i / (N + 1) * Ts[n]
    while (t > Ts[j]) { tj_1 <- Ts[j]; xj_1 <- Xs[j]; yj_1 <- Ys[j]; j <- j + 1; }
    expXs[iS, iA, i] <- (Xs[j] * (t - tj_1) + xj_1 * (Ts[j] - t)) / (Ts[j] - tj_1)
    expYs[iS, iA, i] <- (Ys[j] * (t - tj_1) + yj_1 * (Ts[j] - t)) / (Ts[j] - tj_1)
    i <- i + 1
  }
  
  # The end of loop A
  iA <- iA + 1
}

  # The end of loop S
  iS <- iS + 1
}

# Results

# Checking the distributions of x and y for each s and a by plotting graphs, and their correlations by calculating them
expXYs <- matrix(nrow = length(Ss), ncol = length(As))
iS <- 1
while (iS <= length(Ss)) {
iA <- 1
while (iA <= length(As)) {
  plot(density(expXs[iS, iA, ], bw = "SJ"))
  normalxs <- seq(min(expXs[iS, iA, ]) - 1, max(expXs[iS, iA, ]) + 1, length = 100)
  normalys <- exp(-1 / 2 * normalxs^2) / sqrt(2 * pi)
  lines(normalxs, normalys, col = "red", lwd = 2)
  plot(density(expYs[iS, iA, ], bw = "SJ"))
  normalxs <- seq(min(expYs[iS, iA, ]) - 1, max(expYs[iS, iA, ]) + 1, length = 100)
  normalys <- exp(-1 / 2 / Ss[iS]^2 * normalxs^2) / sqrt(2 * pi) / Ss[iS]
  lines(normalxs, normalys, col = "red", lwd = 2)
  expXYs[iS, iA] <- cor(expXs[iS, iA, ], expYs[iS, iA, ])
  iA <- iA + 1
}
iS <- iS + 1
}

proc_time <- proc.time() - proc_time

# Means and standard deviations of the numbers of SPs and times
expESPNs <- matrix(nrow = length(Ss), ncol = length(As), dimnames = list(rep(0, length(Ss)), rep(0, length(As))))
iS <- 1
while (iS <= length(Ss)) {
  rownames(expESPNs)[iS] <- paste("s = ", Ss[iS], ":", sep = "")
  iA <- 1
  while (iA <= length(As)) {
    if (iS == 1) { colnames(expESPNs)[iA] <- paste("a = ", As[iA], sep = ""); }
    expESPNs[iS, iA] <- paste(mean(expSPNs[iS, iA, ]), "(", sd(expSPNs[iS, iA, ]), ")")
    iA <- iA + 1
  }
  iS <- iS + 1
}

expETs <- matrix(nrow = length(Ss), ncol = length(As), dimnames = list(rep(0, length(Ss)), rep(0, length(As))))
iS <- 1
while (iS <= length(Ss)) {
  rownames(expETs)[iS] <- paste("s = ", Ss[iS], ":", sep = "")
  iA <- 1
  while (iA <= length(As)) {
    if (iS == 1) { colnames(expETs)[iA] <- paste("a = ", As[iA], sep = ""); }
    expETs[iS, iA] <- paste(mean(expTs[iS, iA, ]), "(", sd(expTs[iS, iA, ]), ")")
    iA <- iA + 1
  }
  iS <- iS + 1
}

paste("Time :", proc_time[3], "seconds")
expXYs
expESPNs
expETs
