#install.packages("glmnet")
#install.packages('Rcpp')
#install.packages("latticeExtra")

library('Rcpp')
library("glmnet")
library("mice")
library("Hmisc")
library("dplyr")
library("caret")
library("Matrix")
library("ggplot2")

getwd()
setwd("~/Desktop/NeuroPathPredict/Model_V0/")
df_y<- read.csv("Y_qnp_data_0428.csv")
df_x<- read.csv("X_cov_roi_0515.csv")

rownames(df_x) <- df_x$rois
rownames(df_y) <- df_y$X

df_y<- df_y[-c(1)]
df_x<- df_x[-c(1)]

summary(df_y)


df_x$roi <- as.factor(df_x$roi)
df_x$p.no. <- as.factor(df_x$p.no.)

class(df_x$roi)

df <- cbind(df_y,df_x)
colnames(df)[1] <- "QNP_obs"

y <- as.matrix(df[,1])
x <- as.matrix(df[4:932])
x_train <- model.matrix( ~ ., df[,-1], ignore.intercept = TRUE)

####### Run nested cross-validation with parameter tuning ######

# Define cross-validation folds for outer loop (model assessment)
outer_folds <- createFolds(df$QNP_obs, k = 5)

# Set up grid of hyperparameters to tune
lambda_seq <- 10^seq(-2, 2, length.out = 100)
alpha_seq <- seq(0, 1, length.out = 11)
hyper_grid <- expand.grid(alpha = alpha_seq, lambda = lambda_seq)

# Define cross-validation control for inner loop (hyperparameter tuning)
inner_control <- trainControl(method = "cv", number = 5)

# Define cross-validation control for outer loop (model assessment)
outer_control <- trainControl(method = "cv", index = outer_folds, savePredictions = TRUE)

# Set up nested cross-validation using caret
set.seed(825)
cv_results <- caret::train(x = x_train,
                           y = df$QNP_obs,
                           trControl = outer_control,
                           tuneGrid = hyper_grid,
                           method = "glmnet",
                           family = "gaussian",
                           tuneLength = 10,
                           verboseIter = TRUE,
                           allowParallel = TRUE,
                           trace.it = TRUE)

cv_results
# View best hyperparameters
best_alpha <- cv_results$bestTune$alpha
best_lambda <- cv_results$bestTune$lambda

temp_ncv      = coef(cv_results$finalModel, s = best_lambda)
temp_ncv.data = as.data.frame(summary(temp_ncv))
temp_ncv.name = row.names(temp_ncv)
pred_coef <- temp_ncv.name[temp_ncv.data[-1,1]]

trellis.par.set(caretTheme())
plot(cv_results)

# Train final model with best hyperparameters
final_model <- glmnet(x_train,
                      y = df$QNP_obs,
                      alpha = best_alpha,
                      lambda = best_lambda,
                      family = "gaussian")

print(final_model)
pred_ncv <- predict(final_model, newx = x_train, s = "best_lambda")

###### Assess results of nested cross-validation with factored regions and participants ######

## Plot the predicted values from the nested cross-validation model compared with original values
par(mfrow = c(1, 1))
plot(y, pred_ncv, main = "Nested CV", xlab = "Observed", ylab = "Predicted")
abline(0, 1, col = "red")

## Plot the residuals from the nested cross-validation model compared with original values
par(mfrow = c(1, 1))
plot(y, y - pred_ncv, main = "Nested CV", xlab = "Observed", ylab = "Residuals")
abline(0, 0, col = "red")

##Plot distribution  of observed values and prectided values
par(mfrow = c(1, 1))
hist(y, main = "Observed", xlab = "QNP", col = "blue")
hist(pred_ncv, main = "Nested CV", xlab = "QNP", col = "red", add = TRUE)

##Plot boxplot of observed values and prectided values
par(mfrow = c(1, 2))
boxplot(y, main = "Observed", xlab = "QNP", col = c("blue"))
boxplot(pred_ncv, main = "Predicted", xlab = "QNP", col = c("#ff2a00"))

####### Run nested cross-validation with parameter tuning without factors ######

# Set up nested cross-validation using caret
set.seed(825)
cv_results_a <- caret::train(x = x,
                           y = df$QNP_obs,
                           trControl = outer_control,
                           tuneGrid = hyper_grid,
                           method = "glmnet",
                           family = "gaussian",
                           tuneLength = 10,
                           verboseIter = TRUE,
                           allowParallel = TRUE,
                           trace.it = TRUE)

cv_results_a
# View best hyperparameters
best_alpha_a <- cv_results_a$bestTune$alpha
best_lambda_a <- cv_results_a$bestTune$lambda

temp_ncv_a      = coef(cv_results_a$finalModel, s = best_lambda_a)
temp_ncv_a.data = as.data.frame(summary(temp_ncv_a))
temp_ncv_a.name = row.names(temp_ncv_a)
pred_coef_a <- temp_ncv_a.name[temp_ncv_a.data[-1,1]]

trellis.par.set(caretTheme())
plot(cv_results_a)

# Train final model with best hyperparameters
final_model_a <- glmnet(x,
                      y = df$QNP_obs,
                      alpha = best_alpha_a,
                      lambda = best_lambda_a,
                      family = "gaussian")

print(final_model_a)
pred_ncv_a <- predict(final_model_a, newx = x, s = best_lambda_a)

###### Compare results of nested cross-validation with & without factored regions and participants ######

## Plot the predicted values from the nested cross-validation model compared with original values
par(mfrow = c(2, 1))
plot(y, pred_ncv, main = "Nested CV", xlab = "Observed", ylab = "Predicted")
abline(col = "red")
plot(y, pred_ncv_a, main = "Nested CV_a", xlab = "Observed", ylab = "Predicted")
abline(col = "red")

## Plot the residuals from the nested cross-validation model compared with original values
par(mfrow = c(2, 1))
plot(y, y - pred_ncv, main = "Nested CV", xlab = "Observed", ylab = "Residuals")
abline(h = 0, col = "red")
plot(y, y - pred_ncv_a, main = "Nested CV_a", xlab = "Observed", ylab = "Residuals")
abline(h = 0, col = "red")

##Plot distribution  of observed values and prectided values
par(mfrow = c(1, 1))
hist(y, main = "Observed", xlab = "QNP", col = "blue")
hist(pred_ncv_a, main = "Nested CV", xlab = "QNP", col = "red", add = TRUE)

##Plot boxplot of observed values and prectided values
par(mfrow = c(1, 2))
boxplot(y, main = "Observed", xlab = "QNP", col = c("blue"))
boxplot(pred_ncv_a, main = "Predicted", xlab = "QNP", col = c("#ff2a00"))
