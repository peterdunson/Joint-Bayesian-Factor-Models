library(dplyr)
library(purrr)
library(tidyr)

df_zlen <- content %>%
   select(id, agedays, zlen) %>%
   rename(subj = id, time = agedays, y = zlen) %>%
   filter(!is.na(y))

t_grid <- seq(min(df_zlen$time), max(df_zlen$time), length = 50)

# Fit per-subject smoothing spline and predict on grid
smoothed_list <- df_zlen %>%
   group_by(subj) %>%
   nest() %>%
   mutate(
      fit    = map(data, ~ smooth.spline(.x$time, .x$y, df = 8)),
      y_pred = map(fit, ~ predict(.x, x = t_grid)$y)
   )

# Create all combinations of subj and time
smoothed_long <- expand.grid(
   subj = smoothed_list$subj,
   time = t_grid
) %>%
   arrange(subj, time)

# Flatten the predicted values (row-wise) to fill y
smoothed_long$y <- unlist(smoothed_list$y_pred)

# Relabel subj to consecutive 1...N
smoothed_long <- smoothed_long %>%
   mutate(subj = as.integer(factor(subj)))

# Normalize time to [0,1]
smoothed_long <- smoothed_long %>%
   mutate(time = (time - min(time)) / (max(time) - min(time)))

# Checks
length(unique(smoothed_long$subj)) # 197
max(smoothed_long$subj)            # 197
head(smoothed_long)

