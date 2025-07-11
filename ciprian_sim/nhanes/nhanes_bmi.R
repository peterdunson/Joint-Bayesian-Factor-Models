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

summary(full_mod)




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





dat <- readRDS("nhanes_phthalates_adults_bmi.rds")



set.seed(123)  # for reproducibility

# 1. Load data (already in dat)
# dat <- readRDS("your_data.rds")   # if you need to load

# 2. Center and scale all columns (including the response)
Y <- scale(dat)

# 3. Make a data.frame from scaled data
Y_df <- as.data.frame(Y)

# 4. Define predictors (all except BMI)
predictor_vars <- setdiff(colnames(Y_df), "BMXBMI")

# 5. Train/test split (70% train, 30% test)
n <- nrow(Y_df)
train_idx <- sample(seq_len(n), size = floor(0.7 * n))
train <- Y_df[train_idx, ]
test  <- Y_df[-train_idx, ]

# 6. Fit linear model on train set
full_mod <- lm(BMXBMI ~ ., data = train)

# 7. Predict on test set
pred_bmi <- predict(full_mod, newdata = test)

# 8. Compute test MSE
mse <- mean((test$BMXBMI - pred_bmi)^2)
cat("Test set MSE (standardized BMI):", round(mse, 4), "\n")\







print(setdiff(top5_vars, colnames(Y_df)))




set.seed(125)  # for reproducibility


# Assume dat is already loaded and Y is already centered/scaled:
Y <- scale(dat)
Y_df <- as.data.frame(Y)

# Top 5 predictors (from your screenshot, with stars)
top5_vars <- c("URXECPT", "URXMBP", "URXMHBP", "URXMHP", "URXMZP")

# 1. 70/30 Train/test split
n <- nrow(Y_df)
train_idx <- sample(seq_len(n), size = floor(0.7 * n))
train <- Y_df[train_idx, ]
test  <- Y_df[-train_idx, ]

# 2. Fit linear model using only top 5 predictors
mod_top5 <- lm(BMXBMI ~ ., data = train[, c("BMXBMI", top5_vars)])

summary(mod_top5)

# 3. Predict on test set
pred_bmi <- predict(mod_top5, newdata = test)

# 4. Compute test MSE
mse <- mean((test$BMXBMI - pred_bmi)^2)
cat("Test set MSE (standardized BMI, top 5 predictors):", round(mse, 4), "\n")

