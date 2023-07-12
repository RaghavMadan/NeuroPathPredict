## Elastic net regression to identify significant predictors for NPP pipeline
# Y: quantitative tau pathology calculated over each of the 10 regions
# this model is to check predictors with only the smallest buffer zone radii
# R:4.3.1

library(glmnet) #V:4.1-7
library(Rcpp) #V:1.0.10
library(caret) #V:6.0-94
library(dplyr) #V:1.1.2
library(dotwhisker) #V:0.7.4

## Data prep ##

set.seed(825, sample.kind = "Rounding")
setwd("~/Desktop/NeuroPathPredict/Model_V0/")
df_y <- read.csv("Y_qnp_data_0618.csv")
df_x <- read.csv("X_cov_roi_0618.csv")

# Column check

names(df_x)[1:20]
names(df_y)

dim(df_x)	#7610 930
dim(df_y)	#7610 1

length(unique(df_x[, "roi"])) # 10
length(unique(df_x[, 'p.no.'])) # 761


# prepare group variables for predictors
  df_x$roi   <- as.factor(df_x$roi)
  df_x$p.no. <- as.factor(df_x$p.no.)
  class(df_x$roi) # "factor"
  class(df_x$p.no.) # "factor"
  
# put x and y together
  df <- data.frame(QNP_obs = df_y, df_x)
# check column names:
dim(df) #  7610  931
head(names(df))

# define standardize()

  standardize = function(x) {
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
    list.cols.kp <- list.cols[!grepl("r2|r5|r7|r10|r12.5|r15|r20|r25|r30|r40", list.cols)]
    df.std.bzr1 <- data.frame(QNP_obs = df.std.temp[,1],(df.std[,colnames(df.std) %in% list.cols.kp]))
    dim(df.std.bzr1) # 7610 342
    
# set up training/testing data
    training.samples = defineTrain(nrow(df),0.8)  # get 80% of df as training set # n=7610
    train.data = df.std.bzr1[training.samples==0,]; dim(train.data) # 6072  342
    test.data  = df.std.bzr1[training.samples==1,]; dim(test.data) # 1538  342
    
# check if any row of x has at least one missing value
  sum(is.na(apply(df.std.bzr1[,-c(2,3)], 1, sum)))  # 0
  check.col.na = apply(df.std.bzr1[,-c(2,3)], 2, sum); length(check.col.na); sum(is.na(check.col.na)) # 340 check!   0

# build the model using the training set
	model <- train( QNP_obs ~ ., data = train.data, method = "glmnet",
           	        trControl = trainControl("cv", number = 10),
                    tuneLength = 10, family="gaussian", 
                    VerboseIter = TRUE, allowParallel = TRUE, trace.it = TRUE
	                )

# best tuning parameter
    best.lambda = model$bestTune$lambda # # 0.00352
    best.alpha  = model$bestTune$alpha  # 0.4

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
    # 0.4466465 0.8097785
	  
#	  Identify significant predictors
	  sum(is.na(keep[,1]))	# 0
	  sum(abs(keep[,1])>0.01)	# 737
	  
	  keep.select = keep[ abs(keep[,1]) > 0.0265,1]  # # <<< thr = 0.0265
	  dim(keep.select)	# null
	  length(keep.select)     # 714
	  list1 = names(keep.select)      # this returns significant ROI names as well as some p.no.
	  list2 = keep.select[!grepl("p.no.", list1)]; length(list2) # 14

 # significant ROIs
	  names(list2)
	  # [1] "(Intercept)"           "roiCA4"                "roiDG"                 "bz_CST_R_bz_r1"       
	  # [5] "bz_Yeo2011_7_bz_r1"    "edt_dxFrom_C_FP_L"     "edt_dxFrom_CerebrA_60" "edt_dxFrom_CerebrA_67"
	  # [9] "edt_dxFrom_CerebrA_82" "edt_dxFrom_CerebrA_84" "edt_dxFrom_TR_S_L"     "edt_dxFrom_Yeo2011_3" 
	  # [13] "edt_dxFrom_Yeo2011_6"  "var_Yeo2011_7"

# Subset for lm model
    df.std.lm <- data.frame(QNP_obs = df.std.temp[,1],(df.std[,colnames(df.std) %in% names(list2)]))

# Lm model with only significant predictors
	model.lm.r1 <- lm( QNP_obs ~ ., data = df.std.lm)
	summary(model.lm.r1)
	lm(data.frame(scale(model.lm.r1$model)))

  dwplot(model.lm.r1, dot_args = list(size = 3, pch = 21, fill = "white")) +
    theme_grey(base_size =20) +
    theme(aspect.ratio = 1, legend.position ="none") +
    ylab("Significant predictors") + xlab("Coeffecient values") +
    geom_vline(xintercept = 0, colour = "black", linetype = 2) +
    ggtitle("R1 model")

	  
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
	  