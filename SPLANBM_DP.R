# Code for the article:
# "Semi-parametric Bayesian inference to shape-restricted partially linear
# additive mixed-effects negative binomial models for longitudinal data"
# Authors: Miaojie Xia, Ya Liu, and Jiang Du
# Miaojie Xia and Ya Liu contributed equally.
# Corresponding author: Jiang Du (dujiang@bjut.edu.cn)
# Beijing University of Technology

library(BayesLogit)
library(mvtnorm)

monotone_basis <- function(x, knots) {
  K <- length(knots) - 2
  m <- K + 2
  S <- matrix(1, m, length(x))

  for (j in seq_len(K - 1)) {
    i1 <- x <= knots[j]
    i2 <- x > knots[j] & x <= knots[j + 1]
    i3 <- x > knots[j + 1] & x <= knots[j + 2]
    i4 <- x > knots[j + 2]
    S[j, i1] <- 0
    S[j, i2] <- (x[i2] - knots[j])^2 /
      ((knots[j + 2] - knots[j]) * (knots[j + 1] - knots[j]))
    S[j, i3] <- 1 - (x[i3] - knots[j + 2])^2 /
      ((knots[j + 2] - knots[j + 1]) * (knots[j + 2] - knots[j]))
    S[j, i4] <- 1
  }

  i1 <- x <= knots[K]
  i2 <- x > knots[K] & x <= knots[K + 1]
  i3 <- x > knots[K + 1] & x <= knots[K + 2]
  i4 <- x > knots[K + 2]
  S[K, i1] <- 0
  S[K, i2] <- (x[i2] - knots[K])^2 /
    ((knots[K + 2] - knots[K]) * (knots[K + 1] - knots[K]))
  S[K, i3] <- 1 - (x[i3] - knots[K + 2])^2 /
    ((knots[K + 2] - knots[K + 1]) * (knots[K + 2] - knots[K]))
  S[K, i4] <- 1

  i1 <- x <= knots[2]
  i2 <- x > knots[2]
  S[K + 1, i1] <- 1 - (knots[2] - x[i1])^2 /
    (knots[2] - knots[1])^2
  S[K + 1, i2] <- 1

  i1 <- x <= knots[K + 1]
  i2 <- x > knots[K + 1] & x <= knots[K + 2]
  i3 <- x > knots[K + 2]
  S[K + 2, i1] <- 0
  S[K + 2, i2] <- (x[i2] - knots[K + 1])^2 /
    (knots[K + 2] - knots[K + 1])^2
  S[K + 2, i3] <- 1

  B <- t(S)
  sweep(B, 2, colMeans(B))
}

convex_basis <- function(x, knots) {
  K <- length(knots) - 2
  m <- K + 2
  S <- matrix(1, m, length(x))

  for (j in seq_len(K - 1)) {
    i1 <- x <= knots[j]
    i2 <- x > knots[j] & x <= knots[j + 1]
    i3 <- x > knots[j + 1] & x <= knots[j + 2]
    i4 <- x > knots[j + 2]
    S[j, i1] <- 0
    S[j, i2] <- (x[i2] - knots[j])^3 /
      (3 * (knots[j + 2] - knots[j]) * (knots[j + 1] - knots[j]))
    S[j, i3] <- x[i3] - knots[j + 1] -
      (x[i3] - knots[j + 2])^3 /
      (3 * (knots[j + 2] - knots[j]) * (knots[j + 2] - knots[j + 1])) +
      ((knots[j + 1] - knots[j])^2 -
         (knots[j + 2] - knots[j + 1])^2) /
      (3 * (knots[j + 2] - knots[j]))
    S[j, i4] <- x[i4] - knots[j + 1] +
      ((knots[j + 1] - knots[j])^2 -
         (knots[j + 2] - knots[j + 1])^2) /
      (3 * (knots[j + 2] - knots[j]))
  }

  i1 <- x <= knots[K]
  i2 <- x > knots[K] & x <= knots[K + 1]
  i3 <- x > knots[K + 1]
  S[K, i1] <- 0
  S[K, i2] <- (x[i2] - knots[K])^3 /
    (3 * (knots[K + 2] - knots[K]) * (knots[K + 1] - knots[K]))
  S[K, i3] <- x[i3] - knots[K + 1] -
    (x[i3] - knots[K + 2])^3 /
    (3 * (knots[K + 2] - knots[K]) * (knots[K + 2] - knots[K + 1])) +
    ((knots[K + 1] - knots[K])^2 -
       (knots[K + 2] - knots[K + 1])^2) /
    (3 * (knots[K + 2] - knots[K]))

  i1 <- x <= knots[2]
  i2 <- x > knots[2]
  S[K + 1, i1] <- x[i1] - knots[1] +
    (knots[2] - x[i1])^3 / (3 * (knots[2] - knots[1])^2)
  S[K + 1, i2] <- x[i2] - knots[1]

  i1 <- x <= knots[K + 1]
  i2 <- x > knots[K + 1]
  S[K + 2, i1] <- 0
  S[K + 2, i2] <- (x[i2] - knots[K + 1])^3 /
    (3 * (knots[K + 2] - knots[K + 1])^2)

  B <- t(S)
  H <- cbind(1, x)
  B <- B - H %*% solve(crossprod(H), crossprod(H, B))
  sweep(B, 2, apply(B, 2, function(v) diff(range(v))), "/")
}

shape_basis <- function(Z, shape, K) {
  Z <- as.matrix(Z)
  q <- ncol(Z)
  probs <- seq(0, 1, length.out = K + 2)
  knots <- vapply(seq_len(q), function(j) {
    quantile(Z[, j], probs, type = 1, names = FALSE)
  }, numeric(K + 2))

  basis <- lapply(seq_len(q), function(j) {
    if (shape[j] <= 2) {
      B <- monotone_basis(Z[, j], knots[, j])
    } else {
      B <- convex_basis(Z[, j], knots[, j])
    }
    if (shape[j] %in% c(2, 4)) -B else B
  })

  list(B = do.call(cbind, basis), knots = knots)
}

r_nonnegative_coefficient <- function(a, b, c) {
  x <- max(c, 1e-8)
  power <- a + 1
  for (i in 1:5) {
    u <- runif(1) * exp(-b * (x - c)^2)
    lower <- max(0, c - sqrt(-log(u) / b))
    upper <- c + sqrt(-log(u) / b)
    e <- runif(1)
    x <- (e * upper^power + (1 - e) * lower^power)^(1 / power)
  }
  x
}

stick_weights <- function(v) {
  G <- length(v)
  pi <- numeric(G)
  remaining <- 1
  for (k in seq_len(G - 1)) {
    pi[k] <- v[k] * remaining
    remaining <- remaining * (1 - v[k])
  }
  pi[G] <- remaining
  pi
}

fit_splanbm_dp <- function(y, X, Z, id, shape, K = 3, G = 200,
                           n_iter = 4000, burn = 2000) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  id <- as.integer(factor(id))
  n <- length(y)
  n_subject <- max(id)
  p0 <- ncol(X)
  q <- ncol(Z)
  m <- K + 2

  basis_object <- shape_basis(Z, shape, K)
  B <- basis_object$B
  curved <- shape %in% c(3, 4)
  Z_linear <- sweep(Z[, curved, drop = FALSE], 2,
                    colMeans(Z[, curved, drop = FALSE]))
  W <- cbind(X, Z_linear)
  p <- ncol(W)

  beta_prior_mean <- rep(0, p)
  beta_prior_precision <- diag(0.01, p)
  mu_prior_mean <- 0
  mu_prior_variance <- 100
  sigma_shape <- sigma_rate <- 0.01
  alpha_shape <- alpha_rate <- 0.01
  gamma_shape <- gamma_rate <- 0.01
  r_shape <- r_rate <- 0.01

  beta <- rep(0, p)
  gamma <- rep(0, ncol(B))
  mu <- 0
  sigma2 <- 1
  alpha <- 1
  atom <- rnorm(G)
  v <- c(runif(G - 1), 1)
  pi <- stick_weights(v)
  label <- sample(seq_len(G), n_subject, replace = TRUE)
  epsilon <- atom[label]
  r <- 1

  keep <- n_iter - burn
  beta_chain <- matrix(0, keep, p)
  gamma_chain <- matrix(0, keep, ncol(B))
  epsilon_chain <- matrix(0, keep, n_subject)
  r_chain <- numeric(keep)

  for (iter in seq_len(n_iter)) {
    Wbeta <- drop(W %*% beta)
    Bgamma <- drop(B %*% gamma)
    eta <- Wbeta + Bgamma + epsilon[id]
    kappa <- (y - r) / 2
    omega <- BayesLogit::rpg(n, y + r, eta)

    V_beta <- solve(beta_prior_precision + crossprod(W, W * omega))
    m_beta <- V_beta %*% (
      beta_prior_precision %*% beta_prior_mean +
        crossprod(W, kappa - omega * (Bgamma + epsilon[id]))
    )
    beta <- drop(mvtnorm::rmvnorm(1, m_beta, V_beta))

    Wbeta <- drop(W %*% beta)
    residual <- kappa / omega - Wbeta - epsilon[id] - drop(B %*% gamma)
    for (j in seq_len(ncol(B))) {
      precision <- sum(omega * B[, j]^2)
      score <- sum(omega * residual * B[, j]) + gamma[j] * precision
      old <- gamma[j]
      gamma[j] <- r_nonnegative_coefficient(
        gamma_shape - 1, precision / 2,
        (score - gamma_rate) / precision
      )
      residual <- residual - B[, j] * (gamma[j] - old)
    }
    Bgamma <- drop(B %*% gamma)

    V_mu <- 1 / (G / sigma2 + 1 / mu_prior_variance)
    m_mu <- V_mu * (sum(atom) / sigma2 +
                      mu_prior_mean / mu_prior_variance)
    mu <- rnorm(1, m_mu, sqrt(V_mu))
    sigma2 <- 1 / rgamma(
      1, sigma_shape + G / 2,
      sigma_rate + sum((atom - mu)^2) / 2
    )

    alpha <- rgamma(
      1, alpha_shape + G - 1,
      alpha_rate - sum(log(1 - v[seq_len(G - 1)]))
    )
    for (k in seq_len(G - 1)) {
      v[k] <- rbeta(1, 1 + sum(label == k),
                    alpha + sum(label > k))
    }
    pi <- stick_weights(v)

    tau_atom <- kappa - omega * (Wbeta + Bgamma)
    omega_subject <- as.numeric(rowsum(omega, id))
    tau_subject <- as.numeric(rowsum(tau_atom, id))
    for (k in seq_len(G)) {
      members <- which(label == k)
      V_atom <- 1 / (1 / sigma2 + sum(omega_subject[members]))
      m_atom <- V_atom * (mu / sigma2 + sum(tau_subject[members]))
      atom[k] <- rnorm(1, m_atom, sqrt(V_atom))
    }

    base_eta <- Wbeta + Bgamma
    for (i in seq_len(n_subject)) {
      rows <- which(id == i)
      log_probability <- sapply(seq_len(G), function(k) {
        eta_k <- base_eta[rows] + atom[k]
        log(pi[k]) + sum(dnbinom(
          y[rows], size = r, mu = r * exp(eta_k), log = TRUE
        ))
      })
      probability <- exp(log_probability - max(log_probability))
      label[i] <- sample(seq_len(G), 1, prob = probability)
    }
    epsilon <- atom[label]

    eta <- Wbeta + Bgamma + epsilon[id]
    latent_count <- sapply(y, function(y_i) {
      if (y_i == 0) 0 else {
        prob <- r / (r + seq_len(y_i) - 1)
        sum(runif(y_i) < prob)
      }
    })
    r <- rgamma(
      1, r_shape + sum(latent_count),
      r_rate + sum(log1p(exp(eta)))
    )

    if (iter > burn) {
      s <- iter - burn
      beta_chain[s, ] <- beta
      gamma_chain[s, ] <- gamma
      epsilon_chain[s, ] <- epsilon
      r_chain[s] <- r
    }
  }

  f_chain <- array(0, c(keep, n, q))
  curved_index <- 0
  for (j in seq_len(q)) {
    columns <- (m * (j - 1) + 1):(m * j)
    f_chain[, , j] <- tcrossprod(gamma_chain[, columns], B[, columns])
    if (curved[j]) {
      curved_index <- curved_index + 1
      f_chain[, , j] <- f_chain[, , j] +
        tcrossprod(beta_chain[, p0 + curved_index], Z_linear[, curved_index])
    }
  }

  list(
    beta_chain = beta_chain[, seq_len(p0), drop = FALSE],
    r_chain = r_chain,
    gamma_chain = gamma_chain,
    epsilon_chain = epsilon_chain,
    f_chain = f_chain,
    knots = basis_object$knots
  )
}
