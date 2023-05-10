#install.packages("glmnet")
#install.packages('Rcpp')
#install.packages("mice")

library('Rcpp')
library("glmnet")
library("mice")

getwd()
setwd("~/Desktop/NeuroPathPredict/Model_V0/")
df_y<- read.csv("Y_qnp_data.csv")
df_x<- read.csv("X_cov_roi_mean_val.csv")

rownames(df_x) <- df_x$rois
rownames(df_y) <- df_y$UWA.

df_y<- df_y[-c(1)]
df_x<- df_x[-c(1)]


summary(t(df_y))

# compute each row's mean using mean() function
rowMeans(df_y["MFG",],na.rm=TRUE)
m <- c()
for(i in rownames(df_y)){
  # compute mean for all columns
  mean_value <- rowMeans(df_y[i,],na.rm = TRUE)
  m <- append(m,mean_value)
}

# adding row names to matrix
a <- matrix(m,ncol=1)
rownames(a) <- rownames(df_y)
a

df_y["MFG",][is.na(df_y["MFG",])] <- a["MFG",]
df_y["SMTG",][is.na(df_y["SMTG",])] <- a["SMTG",]
df_y["IPL",][is.na(df_y["IPL",])] <- a["IPL",]
df_y["AM",][is.na(df_y["AM",])] <- a["AM",]
df_y["CA1_L",][is.na(df_y["CA1_L",])] <- a["CA1_L",]
df_y["CA1_R",][is.na(df_y["CA1_R",])] <- a["CA1_R",]
df_y["CA3_L",][is.na(df_y["CA3_L",])] <- a["CA3_L",]
df_y["CA3_R",][is.na(df_y["CA3_R",])] <- a["CA3_R",]
df_y["CA4_L",][is.na(df_y["CA4_L",])] <- a["CA4_L",]
df_y["CA4_R",][is.na(df_y["CA4_R",])] <- a["CA4_R",]
df_y["DG_L",][is.na(df_y["DG_L",])] <- a["DG_L",]
df_y["DG_R",][is.na(df_y["DG_R",])] <- a["DG_R",]
df_y["SB_L",][is.na(df_y["SB_L",])] <- a["SB_L",]
df_y["SB_R",][is.na(df_y["SB_R",])] <- a["SB_R",]
df_y["EC_L",][is.na(df_y["EC_L",])] <- a["EC_L",]
df_y["EC_R",][is.na(df_y["EC_R",])] <- a["EC_R",]

##############################

summary(t(df_y))
summary(df_x)


y <- data.matrix(df_y)
x <- data.matrix(df_x)

####### Run 

cv_model <- cv.glmnet(x, y, family = "mgaussian", nfolds = 8, alpha=1, trace.it = TRUE,standardize.response = TRUE)
print(cv_model)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)
plot(cv_model, xvar = "lambda", label = TRUE, type.coef = "2norm")


best_model = glmnet(x,y, alpha = 1, lambda=best_lambda, family = "mgaussian")
#coef(best_model)
predict(best_model, newx = x[1:3,], s = "lambda.min")

model_se = glmnet(x,y, alpha=1, lambda="lambda.", family = "mgaussian")


temp      = coef(cv_model, s = "lambda.min")
print((temp[[1]]@i))
for (j in length(temp)){
  temp.data[j]<- as.data.frame(summary(temp[j[]]))
}
temp.data = as.data.frametemp[][]  # data.frame of dim 4 x 3
temp.name = row.names(temp)               # vector of x variable names, starting with intercept
temp.name[temp.data[-1,1]]                # "mpg"  "wt"   "qsec" (names of staying x variables, excluding intercept)
