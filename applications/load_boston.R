# load_boston.R

#Boston Housing506 neighborhoods × 13 socioeconomic and housing variables


# install.packages("MASS")
library(MASS)

# Boston: 506 × 14, including the response ‘medv’
data("Boston", package="MASS")

# drop the outcome if you only want predictors
df_boston <- Boston[, setdiff(names(Boston), "medv")]

# quick checks
print(dim(df_boston))
print(summary(df_boston))

# optionally scale
df_boston_scaled <- scale(df_boston)
