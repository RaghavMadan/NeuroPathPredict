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
df_x<- read.csv("X_cov_roi_0505.csv")

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
x <- model.matrix(df[,-1])

x_train <- model.matrix( ~ ., df[,-1], ignore.intercept = TRUE)

class(x[2])
class(y[2])

####### Run nested cross-validation with parameter tuning ######

# Define cross-validation folds for outer loop (model assessment)
outer_folds <- createFolds(df$QNP_obs, k = 10)

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
                           allowParallel = TRUE)

cv_results
# View best hyperparameters
best_alpha <- cv_results$bestTune$alpha
best_lambda <- cv_results$bestTune$lambda

temp_ncv      = coef(cv_results$finalModel, s = best_lambda)
temp_ncv.data = as.data.frame(summary(temp_ncv))
temp_ncv.name = row.names(temp_ncv)           
temp_ncv.name[temp_ncv.data[-1,1]]  

trellis.par.set(caretTheme())
plot(cv_results)
densityplot(cv_results, pch = "|")

# Train final model with best hyperparameters
final_model <- glmnet(x = as.matrix(df[2:978]),
                      y = df$QNP_obs,
                      alpha = best_alpha,
                      lambda = best_lambda,
                      family = "gaussian")

print(final_model)
predict(final_model, newx = x[1:10,], s = "best_lambda")
predict(final_model, newx = x[761:770,], s = "best_lambda")

# Plot the results
ggplot(cv_results) +
  geom_point(aes(x = log(lambda), y = RMSE)) +
  geom_line(aes(x = log(lambda), y = RMSE)) +
  facet_wrap(~ alpha, ncol = 3, scales = "free_x") +
  xlab("Log(lambda)") +
  ylab("RMSE") +
  ggtitle("Nested Cross-Validation Results")

######################End nested cv###########
### original CV model ################

cv_model <- cv.glmnet(x_train,y, family = "gaussian", nfolds = 10, alpha=1, trace.it = TRUE, standardize = FALSE)
print(cv_model)

best_lambda <- cv_model$lambda.min
log(best_lambda)

se_lambda <- cv_model$lambda.1se
se_lambda

plot(cv_model)

best_model = glmnet(x_train,y, alpha = 1, lambda=best_lambda, family = "gaussian")
coef_bm <- as.matrix(coef(best_model))
predict(best_model, newx = x_train[1:10,], s = "lambda.min")
predict(best_model, newx = x[761:770,], s = "lambda.min")

model_se = glmnet(x,y, alpha=0.1, lambda=se_lambda, family = "gaussian")
coef(model_se)
predict(model_se, newx = x[1:5,], s = "lambda.1se")

beta_best_model <- best_model[["beta"]]



temp      = coef(cv_model, s = "lambda.min")
temp.data = as.data.frame(summary(temp)) 
temp.name = row.names(temp)           
temp.name[temp.data[-1,1]]                
