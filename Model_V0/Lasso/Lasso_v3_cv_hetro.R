library('Rcpp')
library("glmnet")
library("mice")
library("Hmisc")
library("dplyr")
library("caret")
library("Matrix")
library("ggplot2")
library("MASS")

getwd()
setwd("~/Desktop/NeuroPathPredict/Model_V0/")

set.seed(825)
df_y<- read.csv("Y_qnp_data_0428.csv")
df_x<- read.csv("X_cov_roi_0515.csv")

rownames(df_x) <- df_x$rois
rownames(df_y) <- df_y$X

df_y<- df_y[-c(1)]
df_x<- df_x[-c(1)]

##log transform df_y##
df_y_ln <- log(df_y)
b <- boxcox(lm(df_y$MFG~1))
l <-b$x[which.max(b$y)]
df_y_bc <- (df_y^l -1)/l

## plot distribution of df_y##
par(mfrow = c(1, 3))
hist(df_y, main = "Observed", xlab = "QNP", col = "blue", breaks = 20)
hist(df_y_ln, main = "Log", xlab = "QNP", col = "blue", breaks = 20)
hist(df_y_bc, main = "Box-Cox", xlab = "QNP", col = "blue", breaks = 20)

df_x$roi <- as.factor(df_x$roi)
df_x$p.no. <- as.factor(df_x$p.no.)

df <- cbind(df_y_ln,df_x)
colnames(df)[1] <- "QNP_obs"

y <- as.matrix(df[,1])
x <- as.matrix(df[4:932])
x_train <- model.matrix( ~ ., df[,-1], ignore.intercept = TRUE)

####### Run nested cross-validation with parameter tuning  for log(y)######

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

plot(cv_results)

# Extract the hyperparameters and performance metrics from cv_results
hyperparams <- data.frame(alpha = cv_results$results$alpha,
                          lambda = cv_results$results$lambda)
performance <- data.frame(RMSE = cv_results$results$RMSE)


# Combine the hyperparameters and performance metrics into a single data frame
data <- cbind(hyperparams, performance)

# Create a scatter plot of the performance metrics against the hyperparameters

ggplot(data, aes(x = log(lambda), y = RMSE, color = factor(alpha))) +
  geom_point() +
  scale_color_discrete(name = "Alpha") +
  labs(x = "Log(Lambda)", y = "RMSE") +
  ggtitle("Hyperparameter Comparison") +
  theme_bw()

# Train final model with best hyperparameters
final_model <- glmnet(x_train,
                      y = df$QNP_obs,
                      alpha = best_alpha,
                      lambda = best_lambda,
                      family = "gaussian")

print(final_model)
pred_ncv <- predict(final_model, newx = x_train, s = best_lambda)

###### Assess results of nested cross-validation with factored regions and participants ######

## Plot the predicted values from the nested cross-validation model compared with original values
par(mfrow = c(1, 1))
plot(y, pred_ncv, main = "Nested CV", xlab = "Observed", ylab = "Predicted")
abline(col = "red")
abline(lm(pred_ncv ~ y), col = "#4c00ff")

## Plot the residuals from the nested cross-validation model compared with original values
par(mfrow = c(1, 1))
plot(y, y - pred_ncv, main = "Nested CV", xlab = "Observed", ylab = "Residuals")
abline(0, 0, col = "red")

##Plot distribution  of observed values and prectided values
par(mfrow = c(1, 1))
hist(y, main = "Observed", xlab = "QNP", col = "blue", breaks = 50)
hist(pred_ncv, main = "Nested CV", xlab = "QNP", col = "red", add = TRUE, breaks =50)

##Plot boxplot of observed values and prectided values
par(mfrow = c(1, 2))
boxplot(y, main = "Observed", xlab = "QNP", col = c("blue"))
boxplot(pred_ncv, main = "Predicted", xlab = "QNP", col = c("#ff2a00"))

####### Run nested cross-validation with parameter tuning  for bc(y)######

df_bc <- cbind(df_y_bc,df_x)
colnames(df_bc)[1] <- "QNP_obs"
y_bc <- df_bc[,1]
x_bc <- model.matrix( ~ ., df_bc[,-1], ignore.intercept = TRUE)

# Set up nested cross-validation using caret
cv_results_bc <- caret::train(x = x_bc,
                           y = df_bc$QNP_obs,
                           trControl = outer_control,
                           tuneGrid = hyper_grid,
                           method = "glmnet",
                           family = "gaussian",
                           tuneLength = 10,
                           verboseIter = TRUE,
                           allowParallel = TRUE,
                           trace.it = TRUE)

cv_results_bc
# View best hyperparameters
best_alpha_bc <- cv_results_bc$bestTune$alpha
best_lambda_bc <- cv_results_bc$bestTune$lambda

temp_ncv_bc      = coef(cv_results_bc$finalModel, s = best_lambda_bc)
temp_ncv_bc.data = as.data.frame(summary(temp_ncv_bc))
temp_ncv_bc.name = row.names(temp_ncv_bc)
pred_coef_bc <- temp_ncv_bc.name[temp_ncv_bc.data[-1,1]]

trellis.par.set(caretTheme())
plot(cv_results_bc)

# Train final model with best hyperparameters
final_model_bc <- glmnet(x_bc,
                      y = df$QNP_obs,
                      alpha = best_alpha_bc,
                      lambda = best_lambda_bc,
                      family = "gaussian")

print(final_model_bc)
pred_ncv_bc <- predict(final_model_bc, newx = x_bc, s = best_lambda_bc)

###### Compare results of nested cross-validation b/w log(y) and bc(y) ######

## Plot the predicted values from the nested cross-validation model compared with original values
par(mfrow = c(2, 1))
plot(y, pred_ncv, main = "Nested CV", xlab = "Observed", ylab = "Predicted")
abline(col = "red")
plot(y_bc, pred_ncv_bc, main = "Nested CV_bc", xlab = "Observed", ylab = "Predicted")
abline(col = "red")

## Plot the residuals from the nested cross-validation model compared with original values
par(mfrow = c(2, 1))
plot(y, y - pred_ncv, main = "Nested CV", xlab = "Observed", ylab = "Residuals")
abline(h = 0, col = "red")
plot(y_bc, y_bc - pred_ncv_bc, main = "Nested CV_bc", xlab = "Observed", ylab = "Residuals")
abline(h = 0, col = "red")

##Plot distribution  of observed values and prectided values
par(mfrow = c(3, 1))
hist(y_bc, main = "Observed", xlab = "QNP", col = "blue", breaks = 50)
hist(pred_ncv, main = "Nested CV", xlab = "QNP", col = "red", add = FALSE, breaks =50)
hist(pred_ncv_bc, main = "Nested CV_bc", xlab = "QNP", col = "green", add = FALSE, breaks =50)

##Plot boxplot of observed values and prectided values
par(mfrow = c(1, 4))
boxplot(y, main = "Observed", xlab = "QNP", col = c("blue"))
boxplot(y_bc, main = "Observed", xlab = "QNP", col = c("blue"))
boxplot(pred_ncv, main = "Predicted log", xlab = "QNP", col = c("#ff2a00"))
boxplot(pred_ncv_bc, main = "Predicted bc", xlab = "QNP", col = c("#00ff33"))


### identify unique coefficients between the two models ###
unique_coef_log <- setdiff(pred_coef, pred_coef_bc)
unique_coef_bc <- setdiff(pred_coef_bc, pred_coef)

### identify better model ###
# The model with the lowest RMSE is the better model
rmse_log <- sqrt(mean((y - pred_ncv)^2))
rmse_bc <- sqrt(mean((y_bc - pred_ncv_bc)^2))

# The model with the lowest MAE is the better model
mae_log <- mean(abs(y - pred_ncv))
mae_bc <- mean(abs(y_bc - pred_ncv_bc))

# The model with the lowest MAPE is the better model
mape_log <- mean(abs((y - pred_ncv)/y))
mape_bc <- mean(abs((y_bc - pred_ncv_bc)/y_bc))

# The model with the lowest R2 is the better model
r2_log <- cor(y, pred_ncv)^2
r2_bc <- cor(y_bc, pred_ncv_bc)^2

##### BC(y) is very slightly better than log(y) #####