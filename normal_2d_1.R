# Initialize alpha
alpha <- runif(1)

# The number of mu samples we want
xyn <- 10000

# Algorithm

# We double xyn to generate twice as many samples as need
# to use the second half later
xyn <- xyn * 2

# 1 Constants
bp = 1 / (1 + alpha)
bm = 1 / (1 - alpha)
c = bp * bm

# 2 Initial values
tk <- replicate(xyn, 0)
xk <- replicate(xyn, 0)
yk <- replicate(xyn, 0)
thetaxk <- replicate(xyn, 0)
thetayk <- replicate(xyn, 0)
k <- 0
t <- 0
x <- 0
y <- 0
thetax <- 1
thetay <- 1

# while loop
while (k < xyn) {

# 3
ax <- c * thetax * (x - alpha * y)
ay <- c * thetay * (y - alpha * x)

# 4
if (thetax == thetay) {
  b <- bp
}
else {
  b <- bm
}

# 5
a = max(0, ax)
ab = a / b
taux <- sqrt(ab * ab - 2 * log(runif(1)) / b) - ab + max(0, -ax / b)
a = max(0, ay)
ab = a / b
tauy <- sqrt(ab * ab - 2 * log(runif(1)) / b) - ab + max(0, -ay / b)

#6
tau = min(taux, tauy)

#7
k <- k + 1
t <- tk[k] <- t + tau
x <- xk[k] <- x + thetax * tau
y <- yk[k] <- y + thetay * tau
if (taux == tau) {
  thetax = -thetax
}
else {
  thetay = -thetay
}
thetaxk[k] <- thetax
thetayk[k] <- thetay

# the end of while loop   
}

# Generate samples
xs <- replicate(xyn, 0)
ys <- replicate(xyn, 0)
j <- 1
i <- 1
t <- tk[1]
while (i <= xyn) {
  t <- t + 1 / (xyn + 1) * (tk[xyn] - tk[1])
  while (tk[j] < t) {
    j <- j + 1
  }
  xs[i] <- xk[j] - thetaxk[j - 1] * (tk[j] - t)
  ys[i] <- yk[j] - thetayk[j - 1] * (tk[j] - t)
  i <- i + 1
}

# Checks

# We ignore the first half of samples
xs <- xs[(xyn / 2 + 1):xyn]
ys <- ys[(xyn / 2 + 1):xyn]

# This graph should be almost linear
# qqnorm shows how well data fits the standard normal distribution
# by plotting data quantiles against theoretical ones
# We normalize xs and ys using their sample means and standard deviations  
qqnorm((xs - mean(xs)) / sd(xs))
qqnorm((ys - mean(ys)) / sd(ys))

# We also report alpha and the sample correlation
c(alpha, cor(xs, ys))
