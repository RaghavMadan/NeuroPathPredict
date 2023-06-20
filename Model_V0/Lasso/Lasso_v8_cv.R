## Lasso regression with cross validation to identify significant predictors for NPP pipeline
# Y: quantitative tau pathology calculated over each of the 10 regions

library("glmnet")
library('Rcpp')
library("mice")
library("Hmisc")
library('caret')

## Data prep ##

set.seed(825)

setwd("~/Desktop/NeuroPathPredict/Model_V0/")
df_y<- read.csv("Y_qnp_data_0618.csv")
df_x<- read.csv("X_cov_roi_0618.csv")

# Column check

names(df_x)[1:20]
names(df_y)

length((df_y[,1]))
length(unique(df_x[,'roi']))
length(unique(df_x[,'p.no.']))

# prepare group variables for predictors
  df_x$roi   <- as.factor(df_x$roi )  ; class(df_x$roi )
  df_x$p.no. <- as.factor(df_x$p.no.) ; class(df_x$p.no.)

# put x and y together
  df <- data.frame(QNP_obs = df_y, df_x)	# put x and y together, excluding their first columns, skipping creating x and y
	# check column names:
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

# set up training/testingi data

    training.samples = defineTrain(nrow(df),0.8)  # get 80% of df as training set
          
    #train.data = df.std[training.samples==0,-1]; dim(train.data)   # train.data has X only, no y	# 6072 x 930 -- approx 80% of data
    #test.data  = df.std[training.samples==1,-1]; dim(test.data)    # test.data  has X only, no y	# 1538 x 930 -- approx 20% of data
                                                                         # n=7610
    train.data = df.std[training.samples==0,]
    test.data  = df.std[training.samples==1,]
    
# check if any row of x has at least one missing value
  sum(is.na(apply(df.std[,-c(1,2)], 1, sum)))  # 0

  check.col.na = apply(df.std[,-c(2,3)], 2, sum); length(check.col.na); sum(is.na(check.col.na)) # 978 check!   0

# build the model using the training set
    #temp.y = df.std[training.samples==0,1]
	model <- train( QNP_obs ~ ., data = train.data, method = "glmnet", 
           	        trControl = trainControl("cv", number = 10),
                    tuneLength = 10, family="gaussian", 
                    VerboseIter = TRUE, allowParallel = TRUE, trace.it = TRUE
	                )

# best tuning parameter
    best.lambda = model$bestTune$lambda
    best.alpha  = model$bestTune$alpha

# coefficient of the final model. You need to specify the best lambda
	  keep = coef(model$finalModel, s = best.lambda)
      keep.data = as.data.frame(summary(keep))
      keep.name = rownames(keep)
      pred_coef <- keep.name[keep.data[-1,1] != 0]  # 0 is the intercept]]
          keep0 = as.array(as.matrix(keep)); rownames(keep0) = NULL

# make predictions on the test data
	  predictions <- model %>% predict(test.data)

# model performance metrics
	  keep1 = data.frame(
		    RMSE    = RMSE( predictions, test.data[, 'QNP_obs']),
		    Rsquare = R2(   predictions, test.data[, 'QNP_obs'])
		  )	
