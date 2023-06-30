## Lasso regression with cross validation to identify significant predictors for NPP pipeline
# Y: quantitative tau pathology calculated over each of the 10 regions
# this model is to check predictors with only the small buffer zone radii, <=r10
# R:4.3.1

library(glmnet) #V:4.1-7
library('Rcpp') #V:1.0.10
#library("mice")
#library("Hmisc")
#library('lattice')
library('caret') #V:6.0-94
library('dplyr') #V:1.1.2

## Data prep ##

set.seed(825, sample.kind = "Rounding")
setwd("~/Desktop/NeuroPathPredict/Model_V0/")
df_y<- read.csv("Y_qnp_data_0618.csv")
df_x<- read.csv("X_cov_roi_0618.csv")

# Column check

names(df_x)[1:20]
names(df_y)

dim(df_x)	#7610 930
dim(df_y)	#7610 1

length(unique(df_x[,'roi']))
length(unique(df_x[,'p.no.']))


# prepare group variables for predictors
  df_x$roi   <- as.factor(df_x$roi )  ; class(df_x$roi )
  df_x$p.no. <- as.factor(df_x$p.no.) ; class(df_x$p.no.)

# put x and y together
  df <- data.frame(QNP_obs = df_y, df_x)	# put x and y together, excluding their first columns, skipping creating x and y
	# check column names:
	dim(df) #  7610  931
	head(names(df))

# define standardize()

  standardize = function(x){
   	   tm = mean(x, na.rm=T)
	   ts = sd(x, na.rm=T)
	   temp = (x-tm)/ts
	   return(temp)
  }
  
# define function to split data into training and test sets
  # n: number of observations
  # p: proportion of observations in training set

  defineTrain = function(n,p){
                        choice1 = runif(n,0,1)
                        choice2 = rep(0, n)
                        choice2[choice1>p] = 1
                        return(choice2)
  }

# standardize x and y
    df.std.temp = apply(df[,-c(2,3)], 2, function(x) standardize(x))
    df.std      = data.frame(QNP_obs = df.std.temp[,1], roi = df[,2], p.no. = df[,3], df.std.temp[,-1])
    dim(df.std.temp); dim(df.std)  # 7610 x  929 , 931

#Subset for only bz_r1
    list.cols <- names(df_x)
    list.cols.kp <- list.cols[!grepl("r12.5|r15|r20|r25|r30|r40", list.cols)]
    df.std.bzr1 <- data.frame(QNP_obs = df.std.temp[,1],(df.std[,colnames(df.std) %in% list.cols.kp]))
    
# set up training/testing data
    training.samples = defineTrain(nrow(df),0.8)  # get 80% of df as training set # n=7610
    train.data = df.std.bzr1[training.samples==0,]; dim(train.data)
    test.data  = df.std.bzr1[training.samples==1,]; dim(test.data)
    
# check if any row of x has at least one missing value
  sum(is.na(apply(df.std.bzr1[,-c(2,3)], 1, sum)))  # 0
  check.col.na = apply(df.std.bzr1[,-c(2,3)], 2, sum); length(check.col.na); sum(is.na(check.col.na)) # 940 check!   0

# build the model using the training set
	model <- train( QNP_obs ~ ., data = train.data, method = "glmnet", 
           	        trControl = trainControl("cv", number = 10),
                    tuneLength = 10, family="gaussian", 
                    VerboseIter = TRUE, allowParallel = TRUE, trace.it = TRUE
	                )

# best tuning parameter
    best.lambda = model$bestTune$lambda # # 0.00859
    best.alpha  = model$bestTune$alpha  # 0.1

# coefficient of the final model. You need to specify the best lambda
	  keep = coef(model$finalModel, s = best.lambda)
      keep.data = as.data.frame(summary(keep))
      keep.name = row.names(keep)
      pred_coef <- keep.name[keep.data[-1,1]]
          keep0 = as.array(as.matrix(keep)); rownames(keep0) = NULL

# make predictions on the test data
	  predictions <- model %>% predict(test.data)

# model performance metrics
	  keep1 = data.frame(
		    RMSE    = RMSE( predictions, test.data[, 'QNP_obs']),
		    Rsquare = R2(   predictions, test.data[, 'QNP_obs'])
		  )	
	  keep1 #       RMSE   Rsquare
    #1 0.4466176 0.8098096
	  
#	  Identify significant predictors
	  sum(is.na(keep[,1]))	# [1] 0
	  sum(abs(keep[,1])>0.01)	# [1] 743 (SC:773)
	  
	  keep.select = keep[ abs(keep[,1]) > 0.03,1]  # <<< need ,1 ow. just numerics
	  dim(keep.select)	# null
	  length(keep.select)     # 743
	  list1 = names(keep.select)      # this returns significant ROI names as well as some p.no.
	  list2 = keep.select[!grepl("p.no.", list1)]; length(list2) # [1] 26 (SC:36)

 # significant ROIs
	  names(list2)
	  
 # Plot the predicted values from the nested cross-validation model compared with original values
	  
	  par(mfrow = c(1, 1))
	  par(cex.lab = 2, cex.axis = 2, cex.main = 2)
	  plot(test.data$QNP_obs, predictions, main = "CV model: test set", xlab = "Observed", ylab = "Predicted")
	  abline(0, 1, col = "red")
	  abline(lm(predictions ~ test.data$QNP_obs), col = "#4c00ff")
	  
# plot predicted values against observed values and group by roi
	  my_palette <- rainbow(length(unique(test.data$roi)))
	  palette(my_palette)
	  par(mfrow = c(1, 1))
	  par(cex.lab = 2, cex.axis = 2, cex.main = 2)
	  plot(test.data$QNP_obs, predictions, main = "QNP values", xlab = "Observed", 
	       ylab = "Predicted", col = test.data$roi, pch =20)
	  legend("bottomright", legend = unique(test.data$roi) , col = my_palette,  pch = 20, cex = 1)
	  abline(col = "red")
	  abline(lm(predictions ~ test.data$QNP_obs), col = "#4c00ff")
	  
# Plot the residuals from the model by roi
	  par(mfrow = c(1, 1))
	  plot(test.data$QNP_obs, test.data$QNP_obs - predictions, main = " residuals in test set", xlab = "Observed",
	       ylab = "Residuals", col = test.data$roi)
	  legend("bottomright", legend = unique(test.data$roi) ,
	         col = my_palette,  pch = 20, cex = 1)
	  abline(0, 0, col = "red")
	  