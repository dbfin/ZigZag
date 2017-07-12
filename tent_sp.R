# The number of sample points we want
xn = 10000

# The number of switching points we generate
spn = xn * 10

# Switching points
xs <- (1 + sqrt(-4 / pi * log(1 - runif(spn)))) * rep(c(1, -1), spn / 2)
ts <- rep(0, spn)
t <- 0
x <- 0
i <- 1
while (i <= spn) {
  t <- ts[i] <- t + abs(xs[i] - x)
  x <- xs[i]
  i <- i + 1
}

# Check whether we indeed have the required distribution
# We draw a graph of the density of a random sample from the process (blue)
# with the theoretical density (red)
xs1 <- replicate(xn, 0)
j <- 1
theta <- 1
i <- 1
while (i <= xn) {
  t <- i / (xn + 1) * ts[spn]
  while (ts[j] < t) {
    j <- j + 1
    theta <- -theta
  }
  xs1[i] <- xs[j] - theta * (ts[j] - t)
  i <- i + 1
}

d1 <- density(xs1, bw = "SJ") # "SJ" is a suggested kernel parameter
plot(d1)
polygon(d1, col = "blue", border = "blue")
fx <- seq(min(xs1 - 1), max(xs1 + 1), length = 100)
fy <- exp(-pi / 4 * (pmax(1, abs(fx)) - 1)^2) / 4
lines(fx, fy, col = "red", lwd = 2)

# Now we take a sample from the generated switching points
xs2 <- xs[sort(sample(1:spn, xn, FALSE))]

# We plot the empirical density of switchng points (blue)
# against the theoretical density of switching points (green),
# and add the density we need (red)
d2 <- density(xs2, bw = "SJ")
plot(d2)
polygon(d2, col = "blue")
spx <- fx
spy <- exp(-pi / 4 * (pmax(1, abs(fx)) - 1)^2) * pi / 4 * (pmax(1, abs(fx)) - 1)
lines(spx, spy, col = "green", lwd = 2)
lines(fx, fy, col = "red", lwd = 2)
