# load_holzinger.R

#Holzinger–Swineford 1939145 students × 9 cognitive-test scores


library(lavaan)

# Holzinger–Swineford 1939: 301 students × 9 tests
data("HolzingerSwineford1939", package="lavaan")

# take only the first school (grades 7–8) for simplicity
df_hs <- subset(HolzingerSwineford1939, school == "Pasteur")[ , paste0("x",1:9)]

# quick checks
print(dim(df_hs))
print(head(df_hs))

# optionally scale
df_hs_scaled <- scale(df_hs)
