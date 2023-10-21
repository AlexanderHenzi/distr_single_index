# required packages
require(isodistrreg)
require(dfoptim)
require(Rcpp)

# functions
sourceCpp("functions/get_crps.cpp")

#' Transform spherical to cartesian coordinates
#' 
#' @param beta vector of length d-1 with first d-2 entries in [0,pi] and last
#'     entry in [0,2pi].
#'     
#' @return 
#' The corresponding cartesian coordinates in R^(d).
spherical_to_cartesian <- function(beta) {
  c(cos(beta), 1) * c(1, cumprod(sin(beta))) 
}

#' Estimate distributional single index model by least squares approach
#' 
#' @param y response variable (vector of size n)
#' @param x covariates (matrix of size n*d)
#' @param m number of grid points in each dimension of spherical coordinates
#'     (see Details).
#' @param optimize at how many points to do numerical optimization to refine the
#'     estimate of the index (see Details).
#' 
#' @details 
#' The optimal index is searched on a grid of the d dimensional unit sphere. The
#' grid points are m uniformly spaced points in [0,pi] for the first d-2 
#' spherical coordinates and 2m uniformly spaced points in [0,2pi] for the
#' (d-1)th coordinate.
#' 
#' If optimize equals 0, then the index is taken as the grid point achieving the
#' minimal value of the target function. If optimize is a positive number k, 
#' then numerical optimization is applied at the k grid points with the smallest
#' error in a second step to find an optimal solution close to the grid point
#' where the error is minimal. If d = 2, this applies the function optimize,
#' and if d > 2 it calls nmkb from the dfoptim package.
#' 
#' @note 
#' The size of the grid for grid search grows fast, it is 2m^(d-1). It the grid
#' gets too large, it is better to choose a moderate m and increase the number
#' in optimize.
#'
#' @return 
#' The estimated index (unit vector of length d), and the output of the 
#' function for optimization (if optimize is positive).
fit_dim <- function(y, x, m, optimize = 1) {
  # compute grid on sphere and the index with alpha on the grid
  d <- ncol(x)
  n <- nrow(x)
  grd <- data.matrix((do.call(
    expand.grid,
    c(
      rep(list(seq(0, pi, length.out = m)), d - 2),
      list(seq(0, 2 * pi, length.out = 2*m))
    )
  )))
  n_y <- rle(sort(y))$lengths
  
  ones <- rep(1, n)
  # compute error for each grid point
  target <- function(beta) {
    alpha <- spherical_to_cartesian(beta)
    index <- round(x %*% alpha, 14) # to avoid negligible numerical differences
    n_x <- aggregate(data.frame(w = ones), by = list(index = index), sum)$w
    fit <- isodistrreg::idr(y = y, X = data.frame(X = index))
    -c(n_x %*% fit$cdf^2 %*% n_y)
  }
  n_grd <- nrow(grd)
  errs <- numeric(n_grd)
  pb <- txtProgressBar(max = n_grd)
  for (i in seq_len(n_grd)) {
    errs[i] <- target(grd[i, ])
    setTxtProgressBar(pb, i)
  }
  close(pb)
  est <- unname(grd[which.min(errs), ])
  opt_result <- NA
  
  # isolate an optimal gridpoint
  if (optimize > 0) {
    ind <- which(errs %in% sort(errs)[seq_len(optimize)])
    err_opt <- min(errs[ind])
    delta <- pi / (m - 1)
    if (d > 2) {
      for (jo in ind) {
        lwr <- pmax(0, grd[jo, ] - delta)
        upr <- c(
          pmin(pi, grd[jo, -(d - 1)] + delta),
          min(2 * pi, grd[jo, d - 1] + delta)
        )
        opt_result <- dfoptim::nmkb(
          fn = target,
          par = (lwr + upr) / 2,
          lower = lwr,
          upper = upr
        )
        if (opt_result$value < err_opt) {
          err_opt <- opt_result$value
          est <- opt_result$par
        }
      }
    } else {
      for (jo in ind) {
        lwr <- pmax(0, grd[jo, ] - delta)
        upr <- c(
          pmin(pi, grd[jo, -(d - 1)] + delta),
          min(2 * pi, grd[jo, d - 1] + delta)
        )
        opt_result <- optimize(
          f = target,
          lower = lwr,
          upper = upr
        )
        if (opt_result$objective < err_opt) {
          err_opt <- opt_result$objective
          est <- opt_result$minimum
        }
      }
    }
  }
  list(alpha = spherical_to_cartesian(est), opt_result = opt_result)
}

#' Estimate monotone single index model
#' 
#' All parameters are exactly the same as for the distributional single index
#' model; see the documentation of the function above.
fit_msim <- function(y, x, m, optimize = 1) {
  # compute grid on sphere and the index with alpha on the grid
  d <- ncol(x)
  grd <- data.matrix((do.call(
    expand.grid,
    c(
      rep(list(seq(0, pi, length.out = m)), d - 2),
      list(seq(0, 2 * pi, length.out = 2*m))
    )
  )))
  
  # compute error for each grid point
  target <- function(beta) {
    alpha <- spherical_to_cartesian(beta)
    index <- c(x %*% alpha)
    yy <- y[order(index, -y)]
    fit <- isoreg(x = index, y = y)$yf
    -mean(fit^2)
  }
  n_grd <- nrow(grd)
  errs <- numeric(n_grd)
  pb <- txtProgressBar(max = n_grd)
  for (i in seq_len(n_grd)) {
    errs[i] <- target(grd[i, ])
    setTxtProgressBar(pb, i)
  }
  close(pb)
  est <- unname(grd[which.min(errs), ])
  opt_result <- NA
  
  # isolate an optimal gridpoint
  if (optimize > 0) {
    ind <- which(errs %in% sort(errs)[seq_len(optimize)])
    err_opt <- min(errs[ind])
    delta <- pi / (m - 1)
    if (d > 2) {
      for (jo in ind) {
        lwr <- pmax(0, grd[jo, ] - delta)
        upr <- c(
          pmin(pi, grd[jo, -(d - 1)] + delta),
          min(2 * pi, grd[jo, d - 1] + delta)
        )
        opt_result <- dfoptim::nmkb(
          fn = target,
          par = (lwr + upr) / 2,
          lower = lwr,
          upper = upr
        )
        if (opt_result$value < err_opt) {
          err_opt <- opt_result$value
          est <- opt_result$par
        }
      }
    } else {
      for (jo in ind) {
        lwr <- pmax(0, grd[jo, ] - delta)
        upr <- c(
          pmin(pi, grd[jo, -(d - 1)] + delta),
          min(2 * pi, grd[jo, d - 1] + delta)
        )
        opt_result <- optimize(
          f = target,
          lower = lwr,
          upper = upr
        )
        if (opt_result$objective < err_opt) {
          err_opt <- opt_result$objective
          est <- opt_result$minimum
        }
      }
    }
  }
  list(alpha = spherical_to_cartesian(est), opt_result = opt_result)
}

#' Estimate distributional single index model with threshold-weighted CRPS
#' 
#' All parameters except for pfun are the same as for the distributional single
#' index model; see the documentation of the function above.
#' 
#' @param pfun a distribution function, taking a single numeric vector as
#'     argument and returning a vector of the same length.
fit_dim_weighted <- function(y, x, m, pfun, optimize = 1) {
  # compute grid on sphere and the index with alpha on the grid
  d <- ncol(x)
  n <- nrow(x)
  grd <- data.matrix((do.call(
    expand.grid,
    c(
      rep(list(seq(0, pi, length.out = m)), d - 2),
      list(seq(0, 2 * pi, length.out = 2*m))
    )
  )))
  sy <- aggregate(list(ind = seq_len(n)), by = list(y_unique = y), list)
  n_y <- lengths(sy$ind)
  p_y <- integer(n)
  p_y[unlist(sy$ind)] <- rep(seq_along(n_y), times = n_y)
  d_y <- diff(pfun(sy$y_unique))

  # compute error for each grid point
  target <- function(beta) {
    alpha <- spherical_to_cartesian(beta)
    index <- round(x %*% alpha, 14) # to avoid negligible numerical differences
    fit <- isodistrreg::idr(y = y, X = data.frame(X = index))
    weighted_crps <- get_crps(
      cdf = fit$cdf,
      d_y = d_y,
      p_y = p_y,
      p_x = fit$indices
    )
    weighted_crps
  }
  n_grd <- nrow(grd)
  errs <- numeric(n_grd)
  pb <- txtProgressBar(max = n_grd)
  for (i in seq_len(n_grd)) {
    errs[i] <- target(grd[i, ])
    setTxtProgressBar(pb, i)
  }
  close(pb)
  est <- unname(grd[which.min(errs), ])
  opt_result <- NA
  
  # isolate an optimal gridpoint
  if (optimize > 0) {
    ind <- which(errs %in% sort(errs)[seq_len(optimize)])
    err_opt <- min(errs[ind])
    delta <- pi / (m - 1)
    if (d > 2) {
      for (jo in ind) {
        lwr <- pmax(0, grd[jo, ] - delta)
        upr <- c(
          pmin(pi, grd[jo, -(d - 1)] + delta),
          min(2 * pi, grd[jo, d - 1] + delta)
        )
        opt_result <- dfoptim::nmkb(
          fn = target,
          par = (lwr + upr) / 2,
          lower = lwr,
          upper = upr
        )
        if (opt_result$value < err_opt) {
          err_opt <- opt_result$value
          est <- opt_result$par
        }
      }
    } else {
      for (jo in ind) {
        lwr <- pmax(0, grd[jo, ] - delta)
        upr <- c(
          pmin(pi, grd[jo, -(d - 1)] + delta),
          min(2 * pi, grd[jo, d - 1] + delta)
        )
        opt_result <- optimize(
          f = target,
          lower = lwr,
          upper = upr
        )
        if (opt_result$objective < err_opt) {
          err_opt <- opt_result$objective
          est <- opt_result$minimum
        }
      }
    }
  }
  list(alpha = spherical_to_cartesian(est), opt_result = opt_result)
}