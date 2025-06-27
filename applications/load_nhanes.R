download.file(
  "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/PHTHTE_J.xpt",
  destfile = "PHTHTE_J.xpt",
  mode = "wb"
)

library(haven)
phthalates <- read_xpt("PHTHTE_J.xpt")


# Download
download.file(
  "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DEMO_J.xpt",
  destfile = "DEMO_J.xpt",
  mode = "wb"
)

# Read it
library(haven)
demo <- read_xpt("DEMO_J.xpt")

str(phthalates)
head(phthalates)
colnames(phthalates)


phthalate_vars <- c(
  "URXCNP", "URXCOP", "URXECP", "URXECPT", "URXHIBP", "URXMBP",
  "URXMC1", "URXMCOH", "URXMEP", "URXMHBP", "URXMHH", "URXMHHT",
  "URXMHNC", "URXMHP", "URXMIB", "URXMNP", "URXMOH", "URXMONP", "URXMZP"
)

phthalates_selected <- phthalates[, c("SEQN", phthalate_vars)]


phthalates_cc <- phthalates_selected[complete.cases(phthalates_selected[, phthalate_vars]), ]



# Merge by SEQN
merged <- merge(phthalates_cc, demo[, c("SEQN", "RIDAGEYR")], by = "SEQN")


phthalates_adults <- merged[merged$RIDAGEYR >= 18, ]


colSums(is.na(phthalates_adults))






dat <- phthalates_adults

# This will open a window and plot the scatterplot matrix (may be slow for 19 variables)
pairs(dat, pch = 20, cex = 0.5, main = "Scatterplot Matrix: NHANES 2017-18 Phthalates")


getwd()

