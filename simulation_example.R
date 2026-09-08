# Simulation example for the SPLANBM-DP model.
# See the article information in SPLANBM_DP.R.

source("SPLANBM_DP.R")

set.seed(2026)

n_subject <- 50
n_time <- 10
n <- n_subject * n_time
id <- rep(seq_len(n_subject), each = n_time)

X <- cbind(x1 = rnorm(n), x2 = rnorm(n))
Z <- cbind(z1 = runif(n), z2 = runif(n))
beta_true <- c(-0.5, 0.5)
r_true <- 2

g1_true <- 1 / (1 + exp(-20 * Z[, 1] + 10)) -
  log((1 + exp(10)) / (1 + exp(-10))) / 20
g2_true <- 3 * (2 * Z[, 2] - 1)^2 - 1

component <- rbinom(n_subject, 1, 0.75)
epsilon_true <- rnorm(n_subject, ifelse(component == 1, 1, -1), 0.5)

eta <- drop(X %*% beta_true) + g1_true + g2_true + epsilon_true[id]
y <- rnbinom(n, size = r_true, mu = r_true * exp(eta))

fit <- fit_splanbm_dp(
  y = y,
  X = X,
  Z = Z,
  id = id,
  shape = c(1, 3),
  K = 3,
  G = 200,
  n_iter = 4000,
  burn = 2000
)

estimates <- c(
  beta1 = mean(fit$beta_chain[, 1]),
  beta2 = mean(fit$beta_chain[, 2]),
  r = mean(fit$r_chain)
)
print(round(estimates, 3))

dev.new(width = 8, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
true_functions <- list(g1_true, g2_true)

for (j in 1:2) {
  order_j <- order(Z[, j])
  posterior_mean <- colMeans(fit$f_chain[, , j])
  lower <- apply(fit$f_chain[, , j], 2, quantile, 0.025)
  upper <- apply(fit$f_chain[, , j], 2, quantile, 0.975)

  plot(Z[order_j, j], posterior_mean[order_j], type = "n",
       xlab = paste0("z", j), ylab = paste0("g", j, "(z)"),
       ylim = range(lower, upper, true_functions[[j]]))
  polygon(
    c(Z[order_j, j], rev(Z[order_j, j])),
    c(lower[order_j], rev(upper[order_j])),
    col = "grey85", border = NA
  )
  lines(Z[order_j, j], posterior_mean[order_j], lwd = 2)
  lines(Z[order_j, j], true_functions[[j]][order_j],
        lwd = 2, lty = 2, col = "red")
  legend("topleft", c("Posterior mean", "True function"),
         lty = c(1, 2), lwd = 2, col = c("black", "red"), bty = "n")
}
