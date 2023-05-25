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

## function to standardize ##
standardize = function(x){
	tm = mean(x, na.rm=T)
	ts = sd(x, na.rm=T)
	temp = (x-tm)/ts
	return(temp)
}

set.seed(825)
df_y<- read.csv("Y_qnp_data_0524.csv")
df_x<- read.csv("X_cov_roi_0524.csv")

rownames(df_y) <- df_y$X

df_y<- df_y[-c(1)]
df_x<- df_x[-c(1)]

df_x$roi <- as.factor(df_x$roi)
df_x$p.no. <- as.factor(df_x$p.no.)
df_x$HC <- as.factor(df_x$HC)

X_s <- df_x[,4:932]
X_std = apply(X_s, 2, function(x) standardize(x))
Y_std = standardize(as.matrix(df_y))
Dat = cbind(Y_std, df_x[,1:3], X_std)
colnames(Dat)[1] <- "QNP_obs"


####### Run nested cross-validation with parameter tuning  for log(y)######

# Define cross-validation folds for outer loop (model assessment)
outer_folds <- createFolds(Dat$QNP_obs, k = 5)

# Set up grid of hyperparameters to tune
lambda_seq <- 10^seq(-2, 2, length.out = 100)
alpha_seq <- seq(0, 1, length.out = 11)
hyper_grid <- expand.grid(alpha = alpha_seq, lambda = lambda_seq)

# Define cross-validation control for inner loop (hyperparameter tuning)
inner_control <- trainControl(method = "cv", number = 5)

# Define cross-validation control for outer loop (model assessment)
outer_control <- trainControl(method = "cv", index = outer_folds, savePredictions = TRUE)

# Set up nested cross-validation using caret
cv_results <- caret::train(x = Dat[,2:933],
                           y = Dat[,1],
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
pred_ncv <- predict(final_model, newx = as.matrix(x_train), s = best_lambda)

###### Assess results of nested cross-validation with factored regions and participants ######

## Plot the predicted values from the nested cross-validation model compared with original values

par(mfrow = c(1, 1))
par(cex.lab = 2, cex.axis = 2, cex.main = 2)
plot(y, pred_ncv, main = "Nested CV", xlab = "Observed", ylab = "Predicted")
abline(0, 1, col = "red")
abline(lm(pred_ncv ~ y), col = "#4c00ff")

## plot predicted values against observed values and group by roi
my_palette <- rainbow(length(unique(df_x$roi)))
palette(my_palette)
par(mfrow = c(1, 1))
par(cex.lab = 2, cex.axis = 2, cex.main = 2)
plot(y, pred_ncv, main = "QNP values", xlab = "Observed", ylab = "Predicted", col = df_x$roi, pch =20)
legend("bottomright", legend = unique(df_x$roi) , col = my_palette,  pch = 20, cex = 2)
abline(col = "red")
abline(lm(pred_ncv ~ y), col = "#4c00ff")

## Plot the residuals from the ncv model by roi
par(mfrow = c(1, 1))
plot(y, y - pred_ncv, main = "Nested CV", xlab = "Observed",
 ylab = "Residuals", col = df_x$roi)
 legend("bottomright", legend = unique(df_x$roi) ,
  col = my_palette,  pch = 20, cex = 1.5)
abline(0, 0, col = "red")

##Plot distribution  of observed values and prectided values
par(mfrow = c(1, 1))
hist(y, main = "Observed", xlab = "QNP", col = df_x$roi, breaks = 50)
hist(pred_ncv, main = "Nested CV", xlab = "QNP", 
col = "red", add = TRUE, breaks =50, )

##Plot boxplot of observed values and prectided values
par(mfrow = c(1, 2))
boxplot(y, main = "Observed", xlab = "QNP", col = c("blue"), ylim = range(c(y, pred_ncv)))
boxplot(pred_ncv, main = "Predicted", xlab = "QNP", col = c("#ff2a00"),ylim = range(c(y, pred_ncv)))

## iterate over 1st column to make 10 columns of 761 values each for 10 rois ##
pred_ncv <- as.data.frame(pred_ncv)

## save predicted values to csv ##
write.csv(pred_ncv, file = "df_lasso_pred_0524.csv")

#save pred_coef to csv
write.csv(pred_coef, file = "pred_coef_0524.csv")
