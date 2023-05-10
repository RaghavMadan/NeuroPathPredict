#install.packages("glmnet")
#install.packages('Rcpp')
install.packages("dplyr")

library('Rcpp')
library("glmnet")
library("mice")
library("Hmisc")
library("dplyr")

getwd()
setwd("~/Desktop/NeuroPathPredict/Model_V0/")
df_y<- read.csv("Y_qnp_data_0428.csv")
df_x<- read.csv("X_cov_roi_0428.csv")

rownames(df_x) <- df_x$rois
rownames(df_y) <- df_y$X

df_y<- df_y[-c(1)]
df_x<- df_x[-c(1)]

summary(df_y)

y <- as.matrix(df_y)
x <- as.matrix(df_x)


####### Run 

cv_model <- cv.glmnet(x,y, family = "gaussian", nfolds = 20, alpha=1, trace.it = TRUE)
print(cv_model)

best_lambda <- cv_model$lambda.min
best_lambda

se_lambda <- cv_model$lambda.1se
se_lambda

plot(cv_model)

best_model = glmnet(x,y, alpha = 1, lambda=best_lambda, family = "gaussian")
coef(best_model)
predict(best_model, newx = x[1:10,], s = "lambda.min")
predict(best_model, newx = x[761:770,], s = "lambda.min")

model_se = glmnet(x,y, alpha=1, lambda=se_lambda, family = "gaussian")
coef(model_se)
predict(model_se, newx = x[1:5,], s = "lambda.1se")

beta_best_model <- best_model[["beta"]]



temp      = coef(cv_model, s = "lambda.min")
temp.data = as.data.frame(summary(temp)) 
temp.name = row.names(temp)           
temp.name[temp.data[-1,1]]                
