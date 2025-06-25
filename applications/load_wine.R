# load_wine.R

#Wine178 wine samples × 13 chemical measurements


library(rattle)      # for wine dataset
data(wine, package="rattle")

# wine: 178 × 13 chemical attributes + Class
df_wine <- wine[, -1]  # drop the Class column if you just want numeric X

# quick checks
print(dim(df_wine))
print(head(df_wine))

# optionally scale
df_wine_scaled <- scale(df_wine)
