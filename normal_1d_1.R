# Initialize mu
mu_true <- rnorm(1)

# Initialize sample
n <- 999999
x <- rnorm(n, mean = mu_true, sd = 1)
mun <- 10000 # The number of mu samples we want

# Algorithm

# We double mun to generate twice as many samples as need
# to use the second half later
mun <- mun * 2

# 1 Constants
s <- sum(x)
b <- n + 1

# 2 Initial values
tk <- replicate(mun, 0)
muk <- replicate(mun, 0)
k <- 0
t <- 0
mu <- 0
theta <- 1

# while loop
while (k < mun) {

# 3
a <- theta * ( (n + 1) * mu - s )
ab <- a / b

# 4
if (a < 0) {
  t <- t - ab
  mu <- mu - ab * theta
  a <- 0
  ab <- 0
}

# 5
tau <- sqrt(ab * ab - 2 * log(runif(1)) / b) - ab

#6
k <- k + 1
tk[k] <- t + tau
muk[k] <- mu + theta * tau
t <- tk[k]
mu <- muk[k]
theta <- -theta

# the end of while loop   
}

# Generate mu samples
mus <- replicate(mun, 0)
j <- 1
theta <- 1
i <- 1
while (i <= mun) {
  t <- i / (mun + 1) * tk[mun]
  while (tk[j] < t) {
    j <- j + 1
    theta <- -theta
  }
  mus[i] <- muk[j] - theta * (tk[j] - t)
  i <- i + 1
}

# Checks

# We ignore the first half of mu samples
mus <- mus[(mun / 2 + 1):mun]

qqnorm((mus - mean(mus)) / sd(mus))

# 1) mu_true is the true value of mu
# 2) mu contains the final value of mu, which should be close to mu_true
# because the posterior variance is small
# 3) the mean of mus is the mean of the obtained sample, so it should be almost mu_true
# 4) the standard deviation of mus should be approximately (n+1)^(-1)
c(mu_true, mu, mean(mus), sd(mus))
