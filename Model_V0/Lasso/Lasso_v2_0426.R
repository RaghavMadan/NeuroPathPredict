#install.packages("glmnet")
#install.packages('Rcpp')
install.packages("Hmisc")

library('Rcpp')
library("glmnet")
library("mice")
library("Hmisc")

getwd()
setwd("~/Desktop/NeuroPathPredict/Model_V0/")
df_y<- read.csv("Y_qnp_data_preproc_0426.csv")
df_x<- read.csv("X_cov_roi_mean_val_0426.csv")

rownames(df_x) <- df_x$rois
rownames(df_y) <- df_y$X

df_y<- df_y[-c(1)]
df_x<- df_x[-c(1)]

summary(t(df_y))

#df_yt <- data.frame(t(df_y))
#y_roi <- df_yt$MFG

#df_xt <- data.frame(t(df_x))
#x_roi <- df_xt$MFG

y <- data.matrix(df_y)
x <- data.matrix(df_x)

####### Run 

cv_model <- cv.glmnet(x, y, family = "gaussian", nfolds = 10, alpha=1, trace.it = TRUE,standardize.response = TRUE)
print(cv_model)

best_lambda <- cv_model$lambda.min
best_lambda

se_lambda <- cv_model$lambda.1se
se_lambda

plot(cv_model)

best_model = glmnet(x,y, alpha = 1, lambda=best_lambda, family = "mgaussian")
#coef(best_model)
predict(best_model, newx = x[1:3,], s = "lambda.min")

model_se = glmnet(x,y, alpha=1, lambda=se_lambda, family = "mgaussian")
predict(model_se, newx = x[1:5,], s = "lambda.1se")

beta_best_model <- best_model[["beta"]]



temp      = coef(cv_model, s = "lambda.min")
print((temp[[1]]@i))
for (j in length(temp)){
  temp.data[j]<- as.data.frame(summary(temp[j[]]))
}
temp.data = as.data.frametemp[][]  # data.frame of dim 4 x 3
temp.name = row.names(temp)               # vector of x variable names, starting with intercept
temp.name[temp.data[-1,1]]                # "mpg"  "wt"   "qsec" (names of staying x variables, excluding intercept)
