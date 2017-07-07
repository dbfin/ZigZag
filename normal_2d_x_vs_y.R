# Initialize alpha
alpha <- 0.99

# The number of mu samples we want
xyn <- 1000

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

# We ignore the first half of samples
xk <- xk[(xyn / 2 + 1):xyn]
yk <- yk[(xyn / 2 + 1):xyn]

# xk vs yk
plot(xk, yk, type="n")
lines(xk, yk, type="l")
