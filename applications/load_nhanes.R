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


getwd()

# Update the path if your R working directory is not your Downloads folder
bodymeasures <- read_xpt("/Users/peterdunson/Desktop/NHANES_bodymeasures_17-18.xpt")

# Check what's inside
str(bodymeasures)
colnames(bodymeasures) # Look for "BMXBMI"



str(phthalates)
head(phthalates)
colnames(phthalates)

colnames(demo)

phthalate_vars <- c(
  "URXCNP", "URXCOP", "URXECP", "URXECPT", "URXHIBP", "URXMBP",
  "URXMC1", "URXMCOH", "URXMEP", "URXMHBP", "URXMHH", "URXMHHT",
  "URXMHNC", "URXMHP", "URXMIB", "URXMNP", "URXMOH", "URXMONP", "URXMZP"
)

phthalates_selected <- phthalates[, c("SEQN", phthalate_vars)]


phthalates_cc <- phthalates_selected[complete.cases(phthalates_selected[, phthalate_vars]), ]



# Merge by SEQN
merged <- merge(phthalates_cc, demo[, c("SEQN", "RIDAGEYR" )], by = "SEQN")


phthalates_adults <- merged[merged$RIDAGEYR >= 18, ]


colSums(is.na(phthalates_adults))






dat <- phthalates_adults[, setdiff(names(phthalates_adults), c("SEQN", "RIDAGEYR"))]


saveRDS(dat, file = "nhanes_phthalates_adults.rds")



# Merge in BMI
phthalates_adults_bmi <- merge(
   phthalates_adults,
   bodymeasures[, c("SEQN", "BMXBMI")],
   by = "SEQN",
   all.x = TRUE
)

# Optional: Only keep adults with non-missing BMI
phthalates_adults_bmi <- phthalates_adults_bmi[!is.na(phthalates_adults_bmi$BMXBMI), ]

# Drop SEQN and RIDAGEYR if you want just the analysis vars + BMI
dat_bmi <- phthalates_adults_bmi[, setdiff(names(phthalates_adults_bmi), c("SEQN", "RIDAGEYR"))]

# Save to RDS
saveRDS(dat_bmi, file = "nhanes_phthalates_adults_bmi.rds")


