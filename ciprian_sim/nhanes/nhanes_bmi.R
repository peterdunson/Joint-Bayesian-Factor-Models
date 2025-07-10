# Load data with BMI included
dat <- readRDS("nhanes_phthalates_adults_bmi.rds")

# Identify predictor variables (exclude BMI itself)
predictor_vars <- setdiff(colnames(dat), "BMXBMI")

# 1. Scatterplots: BMI vs. each predictor
par(mfrow = c(4, 5), mar = c(3, 3, 2, 1))  # Adjust layout as needed (for 19 predictors)
for (v in predictor_vars) {
   plot(dat[[v]], dat$BMXBMI,
        xlab = v, ylab = "BMI", main = v, pch = 20, cex = 0.6)
}
par(mfrow = c(1, 1))

# 2. Forward selection (AIC) for linear regression
# Remove rows with missing data
df <- na.omit(dat[, c("BMXBMI", predictor_vars)])

# Fit full and null models
full_mod <- lm(BMXBMI ~ ., data = df)
null_mod <- lm(BMXBMI ~ 1, data = df)

# Forward selection
forward_mod <- step(null_mod, scope = formula(full_mod), direction = "forward", trace = 1)
summary(forward_mod)




# Compute correlation between BMI and each exposure
cor_bmi_exposures <- sapply(predictor_vars, function(v) cor(dat[[v]], dat$BMXBMI, use = "complete.obs"))

# Print as a tidy table
cor_table <- data.frame(
   Predictor = predictor_vars,
   Correlation_with_BMI = round(cor_bmi_exposures, 3)
)
print(cor_table)


cor_table <- cor_table[order(-abs(cor_table$Correlation_with_BMI)), ]
print(cor_table)
