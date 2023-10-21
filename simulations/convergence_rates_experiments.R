# packages
library(truncnorm)
library(RGeode)

# functions
source("functions/estimation.R")

# parameters
pars <- expand.grid(
  sim = seq_len(100),
  n = 2^(8:13),
  distribution = c("normal", "exponential"),
  phi = list(
    pi / 4,
    pi / 3,
    pi / 2,
    c(pi / 4, pi / 2),
    c(pi / 3, pi / 3),
    c(pi / 2, pi / 4)
  ),
  method = c("empirical", "weighted_uniform", "weighted_adapted")
)
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) + 7200

phi <- pars$phi[[id]]
n <- pars$n[id]
distr <- pars$distr[id]
method <- pars$method[id]

if (distr == "normal") {
  dfun <- function(y, x) pnorm(y, mean = x)
  rfun <- function(n, x) rnorm(n, mean = x)
} else if (distr == "exponential") {
  dfun <- function(y, x) pexp(q = y, rate = 1 / x)
  rfun <- function(n, x) rexp(n = n, rate = 1 / x)
}

alpha <- spherical_to_cartesian(phi)
d <- length(alpha)
n_sim_eval <- 5000

# simulate data
set.seed(id - 7200)
x <- matrix(nrow = n, ncol = d, runif(n * d))
index <- c(x %*% alpha)^3
y <- rfun(n, index)

# estimate distributional index model
if (method == "empirical") {
  fit <- fit_dim(y = y, x = x, m = 20, optimize = 1)
} else if (method == "weighted_uniform") {
  if (distr == "normal") {
    pfun <- function(x) punif(x, -10, 10)
    rtfun <- function(n) runif(n, -10, 10)
  } else {
    pfun <- function(x) punif(x, 0, 50)
    rtfun <- function(n) runif(n, 0, 50)
  }
  fit <- fit_dim_weighted(y = y, x = x, m = 20, optimize = 1, pfun = pfun)
} else if (method == "weighted_adapted") {
  if (distr == "normal") {
    pfun <- function(x) {
      out <- (pnorm(x, sd = sqrt(4)) - pnorm(-4, sd = sqrt(4))) / 
        (pnorm(10, sd = sqrt(4)) - pnorm(-4, sd = sqrt(4)))
      out[x < -4] <- 0
      out[x > 10] <- 1
      out
    }
    rtfun <- function(n) rtruncnorm(n, sd = sqrt(4), a = -10, b = 10)
  } else {
    pfun <- function(x) {
      out <- pgamma(x, shape = 3, rate = 1) / pgamma(10, shape = 3, rate = 1)
      out[x > 10] <- 1
      out
    }
    rtfun <- function(n) rgammatr(n, A = 3, B = 1, range = c(0, 10))
  }
  fit <- fit_dim_weighted(y = y, x = x, m = 20, optimize = 1, pfun = pfun)
}
alpha_hat <- fit$alpha
idr_fit <- idr(
  y = y,
  X = data.frame(index = c(x %*% alpha_hat)^3)
)

# compute errors
x_eval <- matrix(nrow = n_sim_eval, ncol = d, runif(n_sim_eval * d))
index_eval <- c(x_eval %*% alpha)^3
index_hat_eval <- c(x_eval %*% alpha_hat)^3
if (method == "empirical") {
  t_eval <- rfun(n_sim_eval, index_eval)
} else {
  t_eval <- rtfun(n_sim_eval)
}
estimated_cdf_x <- cdf(
  predict(idr_fit, data.frame(index = index_hat_eval), digits = 6),
  thresholds = t_eval
)
true_cdf_x <- outer(
  X = index_eval,
  Y = t_eval,
  FUN = function(X, Y) dfun(Y, X)
)
bundled_error <- sqrt(mean((estimated_cdf_x - true_cdf_x)^2))

z_eval <- runif(n_sim_eval, 0, sum(alpha)^3)
estimated_cdf_z <- cdf(
  predict(idr_fit, data.frame(index = z_eval), digits = 6),
  thresholds = t_eval
)
true_cdf_z <- outer(
  X = z_eval,
  Y = t_eval,
  FUN = function(X, Y) dfun(Y, X)
)
cdf_error <- sqrt(mean((estimated_cdf_z - true_cdf_z)^2))

alpha_error <- sqrt(sum((alpha - alpha_hat)^2))

# export results
out <- pars[id, ]
out$bundled_error <- bundled_error
out$cdf_error <- cdf_error
out$alpha_error <- alpha_error

saveRDS(object = out, file = paste0("simulations_", id, ".rds"))