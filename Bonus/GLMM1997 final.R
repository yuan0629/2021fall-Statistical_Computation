library(ggplot2)
M <- 100
beta <- 5
sigma <- sqrt(0.5)
q <- 10
n <- 15
eps <- 1e-4

# MCEM
N <- 1000
N1 <- 500 #burn in
betafinal1 <- rep(0, M)
sigmafinal1 <- rep(0, M)
maxbeta <- function(beta, tempu) {
  q1 <- q2 <- q3 <- 0
  q1 <- beta * sum(Y * X)
  # for (t in (N1+1):N) { q2 has no relation with beta
  #   q2 <- q2 + matrix(apply(Y, 2, sum), nrow = 1) %*% matrix(tempu[t, ], ncol = 1)
  # }
  # q2 <- q2 / (N - N1)
  for (t in (N1+1):N) {
    q3 <- q3 + sum(log(1 + exp(beta * X + matrix(tempu[t, ], nrow = n, ncol = q, byrow = TRUE))))
  }
  q3 <- q3 / (N - N1)
  #return(q1 + q2 - q3)
  return(q1 - q3)
}
for (m1 in 95:M) {
  set.seed(m1)
  print(m1)
  # data generation
  u <- rnorm(q, 0, sigma)
  X <- p <- Y <- matrix(0, ncol = q, nrow = n)
  for (j in 1:q) {
    X[, j] <- c(1:n) / n
    p[, j] <- exp(beta * X[, j] + u[j]) / (1 + exp(beta * X[, j] + u[j])) 
  }
  for (i in 1:n) {
    for (j in 1:q) {
      Y[i, j] <- rbinom(1, 1, p[i, j])
    }
  }
  m <- 1
  beta0 <- c(2, rep(0, 100))
  sigma0 <- c(0.3, rep(0, 100))
  while(TRUE) {
    # Metropolis
    tempu <- matrix(0, ncol = q, nrow = N)
    for (t in 2:N) {
      for (j in 1:q) {
        uni <- runif(1, 0, 1)
        uj <- rnorm(1, 0, sqrt(sigma0[m]))
        A <- min(1, exp(sum(Y[, j]) * (uj - tempu[t-1, j])) * prod((1+exp(beta0[m]*X[,1]+tempu[t-1, j]))
                                                                   /(1+exp(beta0[m]*X[,1]+uj)))) #xij is same for j
        if (uni < A) {tempu[t, j] <- uj}
        else {tempu[t, j] <- tempu[t-1, j]}
      }
    }
    # Maximize to get beta(search)
    l <- -1000
    for (b in seq(0, 10, by = 0.1)) {
      if (maxbeta(b, tempu) > l) {
        beta0[m+1] <- b
        l <- maxbeta(b, tempu)
      }
    }
    # update to get sigma
    sigma0[m+1] <- 1 / (N - N1) * sum(tempu[((N1+1):N),]^2 / q)
    m <- m + 1
    #print(m)
    if((abs(beta0[m] - beta0[m-1]) < eps) & (abs(sigma0[m] - sigma0[m-1]) < eps)) {break}
  }
  betafinal1[m1] <- beta0[m]
  sigmafinal1[m1] <- sigma0[m]
}


# MCNR
N <- 10000
N1 <- 5000 #burn in
betafinal2 <- rep(0, M)
sigmafinal2 <- rep(0, M)
funcW <- function(beta, tempu) {
  W <- matrix(0, nrow = n, ncol = n)
  mu <- matrix(0, ncol = 1, nrow = n)
  for (i in 1:n) {
    mu1 <- mean(1 / (1 + exp(-beta * X[i,1] - tempu)))
    mu[i, 1] <- mu1
    W[i, i] <- mu1 * (1 - mu1)
  }
  return(list(W = W, mu = mu))
}
for (m1 in 10:M) {
  set.seed(m1)
  print(m1)
  u <- rnorm(q, 0, sigma)
  X <- p <- Y <- matrix(0, ncol = q, nrow = n)
  for (j in 1:q) {
    X[, j] <- c(1:n) / n
    p[, j] <- exp(beta * X[, j] + u[j]) / (1 + exp(beta * X[, j] + u[j])) 
  }
  for (i in 1:n) {
    for (j in 1:q) {
      Y[i, j] <- rbinom(1, 1, p[i, j])
    }
  }
  m <- 1
  beta0 <- c(2, rep(0, 100))
  sigma0 <- c(0.7, rep(0, 100))
  while(TRUE) {
    # Metropolis
    tempu <- matrix(0, ncol = q, nrow = N)
    for (t in 2:N) {
      for (j in 1:q) {
        uni <- runif(1, 0, 1)
        uj <- rnorm(1, 0, sqrt(sigma0[m]))
        A <- min(1, exp(sum(Y[, j]) * (uj - tempu[t-1, j])) * prod((1+exp(beta0[m]*X[,1]+tempu[t-1, j]))
                                                                   /(1+exp(beta0[m]*X[,1]+uj)))) #xij is same for j
        if (uni < A) {tempu[t, j] <- uj}
        else {tempu[t, j] <- tempu[t-1, j]}
      }
    }
    # Newton-Raphson to get beta
    W <- funcW(beta0[m], tempu[((N1+1):N),])
    X1 <- matrix(X[,1], nrow = n)
    beta0[m+1] <- beta0[m] + 1 / (t(X1) %*% W$W %*% X1) * (t(X1) %*% matrix(apply(Y,1,mean) - W$mu, ncol = 1))
    # update sigma
    sigma0[m+1] <- 1 / (N - N1) * sum(tempu[((N1+1):N),]^2 / q)
    m <- m + 1
    #print(m)
    if((abs(beta0[m] - beta0[m-1]) < eps) & (abs(sigma0[m] - sigma0[m-1]) < eps)) {break}
  }
  betafinal2[m1] <- beta0[m]
  sigmafinal2[m1] <- sigma0[m]
}
