# The number of sample points we want
xn = 10000

# The number of switching points we generate
spn = xn * 10

# Switching points
xs <- sqrt(-2 * log(1 - runif(spn))) * rep(c(1, -1), spn / 2)
ts <- rep(0, spn)
t <- 0
x <- 0
i <- 1
while (i <= spn) {
  t <- ts[i] <- t + abs(x) + abs(xs[i])
  x <- xs[i]
  i <- i + 1
}

# Check whether we indeed have the standard normal
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

qqnorm(xs1)

# Now we take switching points
xs2 <- xs[sort(sample(1:spn, xn, FALSE))]

# And we plot their density against the standard normal density
d <- density(xs2)
plot(d)
polygon(d, col = "blue")
normx <- seq(min(xs2), max(xs2), length=100)
lines(normx, dnorm(normx), col="red", lwd=2)
