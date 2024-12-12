#------------- Script for prediction models insomnia symptom trajectories --------------#
#
# Description:  This script does prediction of insomnia symptom trajectories wit both linear and nonlinear models
#
# Authors:      B.Bruinsma
# Date:         December 2024
# Version:      2.0
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                          Settings & Dependencies
#------------------------------------------------------------------------------#
setwd("Folder location")

# Define path of the Rproject to get and save files
save_loc = "Folder location"

# data import
library(haven)
library(readxl)

# data manipulation
library(tidyverse)
library(naniar)
library(reshape2)

# visualizations
library(ggplot2)

# statistics
library(openxlsx)
library(matrixStats)
library(pROC)
library(caret)
library(xgboost)
library(shapr)
library(car)
library(stats)
library(glmnet)

# Other
library(doMC)


#------------------------------------------------------------------------------#
#                              Data Collection
#------------------------------------------------------------------------------#
# # Load the RData file
load("LocationToFile.RData")


#------------------------------------------------------------------------------#
#                       Prep data for prediction models
#------------------------------------------------------------------------------#

# Adapt binary coding: convert 1 to 0 and 2 to 1 in sleep_trajec_binary (0 = resilient sleeper, 1 = insomnia symptoms)
Stat_data$sleep_trajec_binary <- ifelse(Stat_data$sleep_trajec_binary == 1, 0, 1)

# Remove the non-relevant columns of Stat_data for this analysis (eg subject ID)
Stat_data <- Stat_data[, -c(2)]

# Split the data into training and testing sets
set.seed(123) # For reproducibility
indices <- sample(1:nrow(Stat_data), size = 0.8 * nrow(Stat_data))
train_data <- Stat_data[indices, ]
test_data <- Stat_data[-indices, ]

# Convert target variable to "yes" (= poor sleepers) /"no" (= good sleepers)
train_data$sleep_trajec_binary <- factor(ifelse(train_data$sleep_trajec_binary == 1, "yes", "no")
test_data$sleep_trajec_binary <- factor(ifelse(test_data$sleep_trajec_binary == 1, "yes", "no"))

# Verify the levels of the target variable
print(levels(train_data$sleep_trajec_binary))
print(levels(test_data$sleep_trajec_binary))

#------------------------------------------------------------------------------#
#               XGBoost - based on: Habets Biol Psy 2023 analysis
#------------------------------------------------------------------------------#

#Initialize parallel processing
registerDoMC(detectCores()-2)
getDoParWorkers()

## Fit model
hyperparams <- 1000

## set seed list for reproduction
set.seed(101)
seeds <- vector(mode = "list", length = 101)
for(i in 1:100) seeds[[i]] <- sample.int(1000, hyperparams)
seeds[[101]] <- sample.int(1000,1)

adaptControl <- trainControl(method = "adaptive_cv",
                             number = 10, repeats = 10,
                             adaptive = list(min = 10, alpha = 0.05, method = "gls", complete = TRUE), #using adaptive resampling to have more efficient computation, see http://topepo.github.io/caret/adaptive-resampling.html
                             classProbs = TRUE, 
                             summaryFunction = twoClassSummary,
                             search = "random",
                             allowParallel = TRUE,
                             verboseIter = TRUE)


xgbTune <- train(x = train_data[,-1],
                 y = train_data$sleep_trajec_binary,
                 method = "xgbTree",
                 tuneLength = hyperparams,
                 metric = "ROC",
                 trControl = adaptControl, 
                 preProcess = c("zv"),
                 verbose = TRUE)


xgbTune
# Save the model
saveRDS(xgbTune, "filename.RDS")


# IF already done the training, use:
xgbTune <- readRDS("path to file")
xgbTune$finalModel
plot(density(xgbTune$resample$ROC))


# predictions on test set:
test_data$sleep_trajec_binary <- factor(test_data$sleep_trajec_binary, levels = c("no", "yes"))
XGbPred <- predict(xgbTune, test_data[,-1])
cm <- confusionMatrix(XGbPred, test_data$sleep_trajec_binary)

#Get predicted class probabilities so we can build ROC curve.
XGbProbs <- predict(xgbTune, test_data[,-1], type="prob")
head(XGbProbs)

#build ROC curve
XGbROC <- roc(test_data$sleep_trajec_binary, XGbProbs[,"yes"])
auc(XGbROC)

#Variable importance
var_imp <- varImp(xgbTune, scale = FALSE)
var_imp$importance


#plot ROC curve
ggplot(data.frame(specificity = XGbROC$specificities, sensitivity = XGbROC$sensitivities), aes(specificity, sensitivity)) +
  geom_path() +
  scale_x_reverse() +
  geom_abline(intercept = 1, slope = 1, colour='grey') +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  labs(title = paste0("AUROC =", signif(auc(XGbROC), 3)))


#------------------------------------------------------------------------------#
#           XGBoost: SHAP (Shapley Additive Explanation) of variables
#------------------------------------------------------------------------------#

# With SHAP for xgboost package:

dataX <- as.matrix(train_data[,xgbTune$finalModel$feature_names])
shap_values <- shap.values(xgb_model = xgbTune$finalModel, X_train = dataX)
shap_values$mean_shap_score


# Extract variable names and mean SHAP scores
variable_names <- names(shap_values$mean_shap_score)
mean_shap_scores <- shap_values$mean_shap_score

# Combine into a data frame
shap_data <- data.frame(Variable = variable_names, Mean_SHAP_Score = mean_shap_scores)

# Write to Excel
write.xlsx(shap_data, "shap_scores_postDep_281124.xlsx", rowNames = FALSE)

# For the plot filter out the variables that have SHAP = 0 (optional)
shap_plot_filtered <- subset(shap_data_plot, Mean_SHAP_Score > 0)


#------------------------------------------------------------------------------#
#                           Logistic Regression
#------------------------------------------------------------------------------#

# First check multicollinearity among predictors - VIF
#--------------------------------------------------------#

# Fit logistic regression model
model_logit <- glm(
  formula = sleep_trajec_binary ~ .,
  family = binomial(link = "logit"),
  data = Stat_data
)

# Calculate VIF for predictors
vif_values <- car::vif(model_logit)

# Print VIF values for each predictor
cat("Variance Inflation Factor (VIF):\n")
for (i in 1:length(vif_values)) {
  cat(paste(names(vif_values)[i], ": ", vif_values[i], "\n"))
}


#------------------------------------------------------------------------------#
#       Logistic Regression with elastic  - based on Habets 2023
#------------------------------------------------------------------------------#

#Initialize parallel processing
registerDoMC(detectCores()-2)
getDoParWorkers()

hyperparams <- 1000

set.seed(42)
seeds <- vector(mode = "list", length = 101)
for(i in 1:100) seeds[[i]] <- sample.int(1000, hyperparams)
seeds[[101]] <- sample.int(1000,1)

adaptControl <- trainControl(method = "adaptive_cv",
                             number = 10, repeats = 10,
                             adaptive = list(min = 10, alpha = 0.05, method = "gls", complete = TRUE), #using adaptive resampling to have more efficient computation, see http://topepo.github.io/caret/adaptive-resampling.html
                             classProbs = TRUE,
                             summaryFunction = twoClassSummary,
                             preProcOptions = list(), 
                             search = "random",
                             allowParallel = TRUE,
                             verboseIter = TRUE)


elTune <- train(x = train_data[,-1],
                y = train_data[,1],
                method = "glmnet",
                tuneLength = hyperparams,
                metric = "ROC",
                trControl = adaptControl, 
                preProcess = c("nzv", "center", "scale"), # remove near-zero predictors and standardize the predictors
                verbose = TRUE)


saveRDS(elTune, "filename")
elTune <- readRDS("path to file")

## predictions on test set:
elPred <- predict(elTune, test_data[,-1])
cm <- confusionMatrix(elPred, test_data[,1])

## Variable importance (coefficient sorting)
coef <- coef(elTune$finalModel, elTune$bestTune$lambda)
coef_abs <- data.frame(analyte = rownames(coef), coefficient = coef[,1])
coef_abs <- coef_abs %>%
  dplyr::slice(-1) %>% 
  arrange(desc(abs(coefficient)))
coef_abs <- rbind(data.frame(analyte = "intercept", coefficient = coef[1,]), coef_abs) %>% 
  mutate(coefficient = round(coefficient, 3))

# For the plot filter out the variables that have coefficient = 0 (optional)
coef_filtered <- subset(coef_abs, coefficient > 0)

# Write to Excel
write.xlsx(coef_filtered, "Coefficients_LogResEN_postDep_281124.xlsx", rowNames = FALSE)

#Get predicted class probabilities so we can build ROC curve.
elProbs <- predict(elTune, test_data[,-1], type="prob")
head(elProbs)

#build ROC curve
elROC <- roc(test_data[,1], elProbs[,"no"])
auc(elROC)

#plot ROC curve
library(data.table)
ggplot(data.table("1-FPR" = elROC$specificities, TPR = elROC$sensitivities), aes(`1-FPR`, TPR)) +
  geom_path() +
  scale_x_reverse() +
  geom_abline(intercept = 1, slope = 1, colour='grey') +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_classic() +
  labs(title = paste0("AUROC =", signif(auc(elROC), 3)))
