# The convergence of the Zig-Zag algorithm for N((0,0), (1, as, as, s^2)) using the batch method

# The speed of convergence is respresented by the variance of (sum f(X_i) / n - E_pi(f)) * sqrt(n) / sd_pi(f),
# which is estimated using the batch method for discrete MCs

# The values of s to consider
Ss <- c(1, 1)

# The values of a to consider
As <- c(0, 0.25, 0.5)

# The number of batches and also the sample size in each batch for each s and a, so the total sample is n = N * N
N <- 250
n <- N * N

# Functions f to use to evaluate convergence
# There are no relative errors now, so there is no need to consider positive functions, and we use a few simple ones
Fs <- c(
  function (x, y) { x; },
  function (x, y) { y; },
  function (x, y) { x * y; }
)

# E_pi(f)
EPIFs <- c(
  function (s, a) { 0; },
  function (s, a) { 0; },
  function (s, a) { a * s; }
)

# var_pi(f)
VARPIFs <- c(
  function (s, a) { 1; },
  function (s, a) { s^2; },
  function (s, a) { s^2 * (1 + 2 * a^2); }
)

# Statistics we gather
# - one sample for each s and a
expXs <- array(dim = c(length(Ss), length(As), n))
expYs <- array(dim = c(length(Ss), length(As), n))
# - average f(X_i) for each s, f, a and batch
expEFs <- array(rep(0, length(Ss) * length(Fs) * length(As) * N), dim = c(length(Ss), length(Fs), length(As), N))
# - tau for each s, f and a
expTaus <- array(dim = c(length(Ss), length(Fs), length(As)))

proc_time_initial <- proc.time()[3]
proc_time_last <- proc_time_initial
steps_total <- length(Ss) * length(As) * (4 * n)
steps_count <- 0

# Loop S
iS <- 1
while (iS <= length(Ss)) {

# Loop A
iA <- 1
while (iA <= length(As)) {

  # the current s and a
  s <- Ss[iS]
  a <- As[iA]

  # We run the process until we have 4 * n SPs
  # Arrays for the Zig-Zag process
  Xs <- rep(0, 4 * n)
  Ys <- rep(0, 4 * n)
  Ts <- rep(0, 4 * n)

  # Variables
  iSP <- 1 # the index of the next switching point
  theta_x <- 1 # the current theta_x
  theta_y <- 1 # the current theta_y
  x <- 0 # the current x
  y <- 0 # the current y
  t <- 0 # the current t

  # Constants
  ca <- a / s
  cb <- 1 / s^2
  cc <- 1 / (1 - a^2)
  
  # Loop Zig-Zag
  while (iSP <= 4 * n) {

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

    # update t, x and y
    t <- t + tau
    x <- x + theta_x * tau
    y <- y + theta_y * tau

    # collect values
    Xs[iSP] <- x
    Ys[iSP] <- y
    Ts[iSP] <- t

    # report progress
    steps_count <- steps_count + 1
    proc_time <- proc.time()[3]
    if (proc_time - proc_time_last >= 1) {
      print(paste0("s:", iS, "/", length(Ss), " a:", iA, "/", length(As), " SP:", iSP, "/", 4 * n, " ",
                   round(steps_count / steps_total * 100, 1), "% ",
                   round(proc_time - proc_time_initial, 1), " seconds"))
      proc_time_last <- proc_time
    }
    
    # The end of loop Zig-Zag
    iSP <- iSP + 1
    if (tau == tau_x) { theta_x <- -theta_x; }
    else { theta_y <- -theta_y; }
  }

  # Obtain a sample, while also calculating average f(X_i) in batches
  tj_1 <- 0
  xj_1 <- 0
  yj_1 <- 0
  j <- 1
  Tmax <- Ts[length(Ts)]
  i <- 1
  imodN <- 1
  idivN <- 1
  while (i <= n) {
    t <- i / (n + 1) * Tmax
    while (t > Ts[j]) { tj_1 <- Ts[j]; xj_1 <- Xs[j]; yj_1 <- Ys[j]; j <- j + 1; }
    x <- (Xs[j] * (t - tj_1) + xj_1 * (Ts[j] - t)) / (Ts[j] - tj_1)
    expXs[iS, iA, i] <- x
    y <- (Ys[j] * (t - tj_1) + yj_1 * (Ts[j] - t)) / (Ts[j] - tj_1)
    expYs[iS, iA, i] <- y
    iF <- 1
    while (iF <= length(Fs)) {
      expEFs[iS, iF, iA, idivN] <- expEFs[iS, iF, iA, idivN] + Fs[[iF]](x, y) / N
      iF <- iF + 1
    }
    imodN <- imodN + 1
    if (imodN > N) { imodN <- 1; idivN <- idivN + 1; }
    i <- i + 1
  }
  
  # Estimate tau
  iF <- 1
  while (iF <= length(Fs)) {
    expTaus[iS, iF, iA] <- sum((expEFs[iS, iF, iA, ] - EPIFs[[iF]](Ss[iS], As[iA]))^2) / VARPIFs[[iF]](Ss[iS], As[iA])
    iF <- iF + 1
  }

  # The end of loop A
  iA <- iA + 1
}

  # The end of loop S
  iS <- iS + 1
}

# Results

# Final time
proc_time <- proc.time()[3] - proc_time_initial

# Checking the distributions of x and y for each s and a by plotting graphs, and their correlations by calculating them
# We also plot graphs for each f to check how well the obtained tau fits the sample
# So, we plot #Ss * #As * (1 + #Fs) graphs
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
  iF <- 1
  while (iF <= length(Fs)) {
    plot(density(expEFs[iS, iF, iA, ], bw = "SJ"))
    normalxs <- seq(min(expEFs[iS, iF, iA, ]) - 1, max(expEFs[iS, iF, iA, ]) + 1, length = 100)
    v <- expTaus[iS, iF, iA] * VARPIFs[[iF]](Ss[iS], As[iA]) / N
    normalys <- exp(-1 / 2 * (normalxs - EPIFs[[iF]](Ss[iS], As[iA]))^2 / v) / sqrt(2 * pi * v)
    lines(normalxs, normalys, col = "red", lwd = 2)
    iF <- iF + 1
  }
  expXYs[iS, iA] <- cor(expXs[iS, iA, ], expYs[iS, iA, ])
  iA <- iA + 1
}
iS <- iS + 1
}

# For each s we present a table containing functions in rows and the values of a in columns
# In each cell we list the estimation of tau
expFA <- matrix(nrow = length(Ss) * (2 + length(Fs)) - 1, ncol = 1 + length(As))
row <- 1
iS <- 1
while (iS <= length(Ss)) {
  if (iS > 1) {
    expFA[row, ] <- rep("", 1 + length(As))
    row <- row + 1;
  }
  iF <- 1
  while (iF <= length(Fs)) {
    if (iF == 1) {
      expFA[row, 1] <- paste("s = ", Ss[iS], sep = "")
      iA <- 1
      while (iA <= length(As)) {
        expFA[row, 1 + iA] <- paste("a = ", As[iA], sep = "")
        iA <- iA + 1
      }
      row <- row + 1
    }
    expFA[row, 1] <- paste("f #", iF, sep = "")
    iA <- 1
    while (iA <= length(As)) {
      expFA[row, 1 + iA] <- expTaus[iS, iF, iA]
      iA <- iA + 1
    }
    row <- row + 1
    iF <- iF + 1
  }
  iS <- iS + 1
}

paste("Time :", round(proc_time, 2), "seconds")
expXYs
expFA

