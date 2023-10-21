# packages
library(ggplot2)
library(tidyr)
library(dplyr)

# ggplot2 settings, colors
theme_set(theme_bw(base_size = 12))
colpal <- c(
  "#999999",
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#0072B2",
  "#D55E00",
  "#CC79A7",
  "#F0E442"
)


# functions
source("functions/estimation.R")

# data
data <- readRDS("data/house_price_data.rds")
y <- data$price
x <- data.matrix(data[, colnames(data) != "price"])
n <- nrow(data)

# fit distributional single index model
dim <- fit_dim(y, x, m = 120, optimize = 2000)
idr_fit <- idr(y = y, X = data.frame(index = x %*% dim$alpha))
index_dim <- c(x %*% dim$alpha)
n_grid <- 1001
grid_index_dim <- seq(min(index_dim), max(index_dim), length.out = n_grid)
cond_quants <- qpred(
  predict(idr_fit, data = data.frame(index = grid_index_dim)),
  quantiles = seq(0.1, 0.9, 0.2)
)

dim_df_points <- data.frame(
  index = index_dim,
  price = y,
  method = "Distributional single index model (empirical)"
)
dim_df_quants <- data.frame(
  index = grid_index_dim,
  quantile = c(cond_quants),
  tau = factor(rep(seq(0.1, 0.9, 0.2), each = n_grid)),
  method = "Distributional single index model (empirical)"
)

# fit distributional single index model with weighted crps
dim_weighted <- fit_dim_weighted(
  y,
  x,
  m = 120,
  optimize = 2000,
  pfun = identity
)
idr_fit_weighted <- idr(
  y = y,
  X = data.frame(index = x %*% dim_weighted$alpha)
)
index_dim_weighted <- c(x %*% dim_weighted$alpha)
grid_index_dim <- seq(
  min(index_dim_weighted),
  max(index_dim_weighted),
  length.out = n_grid
)
cond_quants <- qpred(
  predict(idr_fit_weighted, data = data.frame(index = grid_index_dim)),
  quantiles = seq(0.1, 0.9, 0.2)
)

dim_weighted_df_points <- data.frame(
  index = index_dim_weighted,
  price = y,
  method = "Distributional single index model (uniform measure)"
)
dim_weighted_df_quants <- data.frame(
  index = grid_index_dim,
  quantile = c(cond_quants),
  tau = factor(rep(seq(0.1, 0.9, 0.2), each = n_grid)),
  method = "Distributional single index model (uniform measure)"
)

# fit single index model for the mean
msim <- fit_msim(y, x, m = 120, optimize = 2000)
index_msim <- c(x %*% msim$alpha)
mono_fit <- isoreg(x = sort(index_msim), y = y[order(index_msim, -y)])

msim_df_points <- data.frame(
  index = mono_fit$x,
  price = mono_fit$y,
  method = "Monotone single index model"
)

save(list = ls(), file = "house_price_application_tmp.rda")

## code the fit as quantile with tau = 0.5 (for simplifying the plotting; the
## curve is an estimate of the conditional mean and not of the median!)
msim_df_mean <- data.frame(
  index = mono_fit$x,
  quantile = mono_fit$yf,
  tau = 0.5,
  method = "Monotone single index model"
)

# get data for quantile model
quantile_data <- readRDS("single_index_quantile_data.rds")
quantile_df_points <- data.frame(
  index = quantile_data[, 7],
  price = quantile_data[, 6] + mean(y),
  method = "Non-crossing quantile regression"
)
quantile_df_quants <- data.frame(
  index = quantile_data[, 7],
  quantile = c(quantile_data[, -c(6, 7)]) + mean(y),
  tau = factor(rep(seq(0.1, 0.9, 0.2), each = n)),
  method = "Non-crossing quantile regression"
)

# combine data, plot
points_df <- rbind(
  dim_df_points,
  dim_weighted_df_points,
  msim_df_points,
  quantile_df_points
)
quants_df <- rbind(
  dim_df_quants,
  dim_weighted_df_quants,
  msim_df_mean,
  quantile_df_quants
)

quantiles_plot <- ggplot() +
  geom_point(
    data = points_df,
    aes(x = index, y = price)
  ) + 
  geom_line(
    data = quants_df,
    aes(x = index, y = quantile, color = tau, group = tau),
    lwd = 0.9
  ) +
  scale_color_manual(values = colpal[1:5]) +
  facet_wrap(.~method, scales = "free_x") +
  theme(legend.position = "none") +
  labs(x = "Index", y = "Price")

pdf(file = "house_price_plot.pdf", width = 8, height = 4)
print(quantiles_plot)
dev.off()

# numbers in the text
ncqr9 <- c(1, -0.395, 6.489, -0.026)
ncqr9 <- ncqr9 / sqrt(sum(ncqr9^2))

alphas <- rbind(
  dim$alpha,
  dim_weighted$alpha,
  msim$alpha,
  ncqr9
)
write.table(
  cbind(
    paste0(c("dim", "dim_uniform", "msim", "ncqr"), " & "),
    matrix(sprintf(alphas, fmt = "%0.3f"), nrow = 4)
  ),
  eol = "\\\\\n",
  sep = " & ",
  row.names = FALSE
)

cors <- cor(
  cbind(quantile_data[, 7], index_dim, index_dim_weighted, index_msim),
  method = "spearman"
)
cors
round(cors, 2)
