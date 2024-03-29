# packages
library(tidyverse)

# load data
fl <- list.files("data/simulation_results/", pattern = ".rds")
data <- vector("list", 21600)
for (j in seq_along(fl)) {
  data[[j]] <- readRDS(paste0("data/simulation_results/", fl[j]))
  if (ncol(data[[j]]) < 8) {
    data[[j]]$method <- "empirical"
  } else {
    data[[j]]$method <- as.character(data[[j]]$method)
  }
}
data <- do.call(rbind, data)

# transform the variable phi
data$phi <- unlist(lapply(
  data$phi, function(x) paste0(round(x, 2), collapse = ", ")
))

# table in the paper
data %>%
  gather(key = "error_measure", value = "error", contains("error")) %>%
  group_by(n, distribution, phi, error_measure, method) %>%
  summarise(error = mean(error)) %>%
  mutate(n = log2(n), error = log2(error)) %>%
  group_by(distribution, phi, error_measure, method) %>%
  nest() %>%
  mutate(
    rate = map_dbl(.x = data, .f = ~coef(lm(.x$error ~ .x$n))[2])
  ) %>%
  select(-data) %>%
  spread(key = "phi", value = "rate") %>%
  select(method, distribution, error_measure, `0.79`, `1.05`, `1.57`, contains(", ")) %>%
  ungroup() %>%
  arrange(method, distribution, error_measure) %>%
  mutate(
    distribution = ifelse(seq_along(distribution) %% 3 == 2, as.character(distribution), ""),
    method = ifelse(seq_along(method) %% 6 == 1, as.character(method), "")
  ) %>%
  mutate_if(is.numeric, function(x) sprintf(x, fmt = "%0.3f")) %>%
  write.table(sep = " & ", eol = " \\\\\n", row.names = FALSE)

# revised table
revised_table <- data %>%
  gather(key = "error_measure", value = "error", contains("error")) %>%
  mutate(n = -log2(n), error = log2(error)) %>%
  group_by(distribution, phi, error_measure, method) %>%
  nest() %>%
  mutate(
    rate = map_dbl(.x = data, .f = ~coef(lm(.x$error ~ .x$n))[2]),
    stderr =
      map_dbl(.x = data, .f = ~summary(lm(.x$error ~ .x$n))$coefficients[2, 2])
  ) %>% 
  select(-data)

revised_table %>%
  mutate(is_empirical = (method == "empirical")) %>%
  group_by(error_measure, is_empirical) %>%
  summarise(
    rate = mean(rate),
    stderr = sqrt(sum(stderr^2)) / sqrt(length(stderr))
  ) %>%
  mutate(distribution = "aggr", error_type = rate) %>%
  arrange(is_empirical)

revised_table %>%
  # bind_rows(aggr) %>%
  mutate(
    ratechr = sprintf(rate, fmt = "%0.2f"),
    stderrchr = sprintf(ceiling(stderr * 100) / 100, fmt = "%0.2f"),
    rate = paste0(ratechr, " (", stderrchr, ")"),
    method = factor(
      method,
      levels = c("empirical", "weighted_adapted", "weighted_uniform"),
      labels = c("Empirical", "Truncated", "Uniform")
    ),
    distribution = factor(
      distribution,
      levels = c("exponential", "normal"),
      labels = c("Exponential", "Gaussian")
    ),
    error_measure = factor(
      error_measure,
      levels = c("alpha_error", "cdf_error", "bundled_error"),
      labels = c("Index", "CDF", "Bundled")
    )
  ) %>%
  select(-ratechr, -stderrchr, -stderr) %>%
  spread(key = "phi", value = "rate") %>%
  select(method, distribution, error_measure, `0.79`, `1.05`, `1.57`, contains(", ")) %>%
  ungroup() %>%
  arrange(method, distribution, error_measure) %>%
  mutate(
    distribution = ifelse(seq_along(distribution) %% 3 == 2, as.character(distribution), ""),
    method = ifelse(seq_along(method) %% 6 == 1, as.character(method), "")
  ) %>%
  write.table(sep = " & ", eol = " \\\\\n", row.names = FALSE)

# average rates over all settings
data %>%
  gather(key = "error_measure", value = "error", contains("error")) %>%
  group_by(n, distribution, phi, error_measure, method) %>%
  summarise(error = mean(error)) %>%
  mutate(n = log2(n), error = log2(error)) %>%
  group_by(distribution, phi, error_measure) %>%
  nest() %>%
  mutate(rate = map_dbl(.x = data, .f = ~coef(lm(.x$error ~ .x$n))[2])) %>%
  group_by(error_measure) %>%
  summarise(rate = mean(rate))

