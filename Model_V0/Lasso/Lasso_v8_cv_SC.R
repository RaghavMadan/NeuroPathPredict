## Lasso regression with cross validation to identify significant predictors for NPP pipeline
# Y: quantitative tau pathology calculated over each of the 10 regions
# R: 4.3.1

library("glmnet")# glmnet_4.1-7
library('Rcpp')  # Rcpp_1.0.10
#library("mice")
library("Hmisc") # Hmisc_5.1-0 
library('caret') # caret_6.0-94
library(dplyr)	 # dplyr_1.1.2

## Data prep ##

set.seed(825)

#setwd("~/Desktop/NeuroPathPredict/Model_V0/")
df_y<- read.csv("Y_qnp_data_0618.csv")
df_x<- read.csv("X_cov_roi_0618.csv")

# Column check

names(df_x)[1:20]
names(df_y)
	
dim(df_x)	#7610 930
dim(df_y)	#7610 1

length(unique(df_x[,'roi']))	#  10
length(unique(df_x[,'p.no.']))	# 761 

# prepare group variables for predictors
  df_x$roi   <- as.factor(df_x$roi )  ; class(df_x$roi )
  df_x$p.no. <- as.factor(df_x$p.no.) ; class(df_x$p.no.)

# put x and y together
  df <- data.frame(QNP_obs = df_y, df_x)	# put x and y together, excluding their first columns, skipping creating x and y
	# check column names:
	head(names(df))
	dim(df)			#  7610  931


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

# standardize x and y, excluding factors (cols 2, 3)
    df.std.temp = apply(df[,-c(2,3)], 2, function(x) standardize(x))
    df.std      = data.frame(QNP_obs = df.std.temp[,1], roi = df[,2], p.no. = df[,3], df.std.temp[,-1])
    dim(df.std.temp); dim(df.std)  # 7610 x  929 , 931

# set up training/testingi data
   set.seed(825)
    training.samples = defineTrain(nrow(df.std),0.8)  # get 80% of df as training set
          
    train.data = df.std[training.samples==0,]; dim(train.data)	# 6154  931
    test.data  = df.std[training.samples==1,]; dim(test.data)	# 1456  931
    
# check if any row of x has at least one missing value
  sum(is.na(apply(df.std[,-c(2,3)], 1, sum)))  # 0   <<<<<<<<<<<<<<<<<<<< has to be 2,3  not 1,2

  check.col.na = apply(df.std[,-c(2,3)], 2, sum); length(check.col.na); sum(is.na(check.col.na)) # 929 check!   0

# build the model using the training set
    #temp.y = df.std[training.samples==0,1]
	model <- train( QNP_obs ~ ., data = train.data, method = "glmnet", 
           	        trControl = trainControl("cv", number = 10),
                    tuneLength = 10, family="gaussian", 
                    VerboseIter = TRUE, allowParallel = TRUE, trace.it = TRUE
	                )

# best tuning parameter
    best.lambda = model$bestTune$lambda # 0.008433909
    best.alpha  = model$bestTune$alpha  # 0.1

# coefficient of the final model. You need to specify the best lambda
	  keep = coef(model$finalModel, s = best.lambda)
      keep.data = as.data.frame(summary(keep))
      keep.name = rownames(keep)
      pred_coef <- keep.name[keep.data[-1,1] != 0]  # 0 is the intercept]]
          keep0 = as.array(as.matrix(keep)); rownames(keep0) = NULL

# make predictions on the test data
	  predictions <- model %>% predict(test.data)  # <<<<<<<<<<<<<< need library(dplyr)
          # predictions <- model$predict(test.data)

# model performance metrics
	  keep1 = data.frame(
		    RMSE    = RMSE( predictions, test.data[, 'QNP_obs']),
		    Rsquare = R2(   predictions, test.data[, 'QNP_obs'])
		  )	

# <<<<<<<<<<<<<<<<<<<<<<<<<<
sum(is.na(keep[,1]))	# [1] 0
sum(abs(keep[,1])>0.01)	# [1] 773

keep.select = keep[ abs(keep[,1]) > 0.15,1]  # <<< need ,1 ow. just numerics
dim(keep.select)	# null
length(keep.select)     # 773
names(keep.select)      # this returns significant ROI names


# after this, you can do...
[1] run another round of elastic net using selected ROIs (ignoring p.no)
[2] run factor analysis
[3] change the threshold in  keep.select = keep[ abs(keep[,1]) > 0.15,1]   # 0.15 -> any number

THEN, go to LUR or UK


