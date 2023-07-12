## Elastic net regression to identify significant predictors for NPP pipeline
# Y: quantitative tau pathology calculated over each of the 10 regions
# R:4.3.1
# Base model

library(glmnet) #V:4.1-7
library(Rcpp) #V:1.0.10
library(caret) #V:6.0-94
library(dplyr) #V:1.1.2
library(dotwhisker) #V:0.7.4

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

length(unique(df_x[,'roi'])) #10
length(unique(df_x[,'p.no.'])) #761

# prepare group variables for predictors
  df_x$roi   <- as.factor(df_x$roi )  ; class(df_x$roi ) #"factor"
  df_x$p.no. <- as.factor(df_x$p.no.) ; class(df_x$p.no.) #"factor"

# put x and y together
  df <- data.frame(QNP_obs = df_y, df_x)
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
    dim(df.std.temp); dim(df.std)  # 7610 x  929 , 7610  x 931

# set up training/testing data
    training.samples = defineTrain(nrow(df),0.8)  # get 80% of df as training set # n=7610
    train.data = df.std[training.samples==0,]; dim(train.data) # 6072  931
    test.data  = df.std[training.samples==1,]; dim(test.data) # 1538  931
    
# check if any row of x has at least one missing value
  sum(is.na(apply(df.std[,-c(2,3)], 1, sum)))  # 0
  check.col.na = apply(df.std[,-c(2,3)], 2, sum); length(check.col.na); sum(is.na(check.col.na)) # 929 check!   0

# build the model using the training set
	model <- train( QNP_obs ~ ., data = train.data, method = "glmnet", 
           	        trControl = trainControl("cv", number = 10),
                    tuneLength = 10, family="gaussian", 
                    VerboseIter = TRUE, allowParallel = TRUE, trace.it = TRUE
	                )

# best tuning parameter
    best.lambda = model$bestTune$lambda # # 0.00359
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
	  keep1
    #       RMSE   Rsquare
    #1 0.4466176 0.8098096
	  
#	  Identify significant predictors
	  sum(is.na(keep[,1]))	# 0
	  sum(abs(keep[,1])>0.01)	# 743
	  
	  keep.select = keep[ abs(keep[,1]) > 0.01,1]  # <<< thr = 0.01
	  dim(keep.select)	# null
	  length(keep.select)     # 743
	  list1 = names(keep.select)      # this returns significant ROI names as well as some p.no.
	  list2 = keep.select[!grepl("p.no.", list1)]; length(list2) # 26

 # significant ROIs
	  names(list2)
	  #capture.output(names(list2), file = "sig_cov_1.txt")
	  #[1] "roiCA4"                   "bz_Buckner2011_3_bz_r10"  "bz_CerebrA_18_bz_r5"     
	  #[4] "bz_CerebrA_5_bz_r1"       "bz_MdLF_R_bz_r10"         "bz_Yeo2011_12_bz_r40"    
	  #[7] "bz_Yeo2011_9_bz_r7"       "bz_Yeo2011_9_bz_r10"      "bz_Yeo2011_9_bz_r12.5"   
	  #[10] "bz_Yeo2011_9_bz_r20"      "bz_Yeo2011_9_bz_r25"      "bz_Yeo2011_9_bz_r30"     
	  #[13] "bz_Yeo2011_9_bz_r40"      "bz_vesRad_thr_0.7_bz_r30" "bz_vesRad_thr_1.0_bz_r40"
	  #[16] "edt_dxFrom_CNIII_R"       "edt_dxFrom_CerebrA_15"    "edt_dxFrom_CerebrA_25"   
	  #[19] "edt_dxFrom_CerebrA_36"    "edt_dxFrom_CerebrA_3"     "edt_dxFrom_CerebrA_60"   
	  #[22] "edt_dxFrom_CerebrA_7"     "edt_dxFrom_CerebrA_82"    "edt_dxFrom_Yeo2011_3"    
	  #[25] "var_CerebrA_5"            "var_F_R"    
	  
######## Second iteration for significant predictor extraction ###########
	
	##Manual removal of bz predictors with larger radius if more than one is picked
	sig.cov.rm <- c("roiCA4","bz_Yeo2011_9_bz_r10",
	                "bz_Yeo2011_9_bz_r12.5", "bz_Yeo2011_9_bz_r15", "bz_Yeo2011_9_bz_r20",
	                "bz_Yeo2011_9_bz_r25", "bz_Yeo2011_9_bz_r30","bz_Yeo2011_9_bz_r40")
	sig.cov1 = names(list2)[!(names(list2) %in% sig.cov.rm)]
	
  length(sig.cov1) #19
  sig.cov1
  df_x.next <- df.std[,colnames(df.std) %in% sig.cov1] 
  df.next <- data.frame(QNP_obs = df.std$QNP_obs, df_x.next)
	  
  # set up training/testing data
  training.samples.next = defineTrain(nrow(df),0.8)  # get 80% of df as training set # n=7610
  train.data.next = df.next[training.samples==0,]; dim(train.data.next) # 6072   20
  test.data.next  = df.next[training.samples==1,]; dim(test.data.next) # 1538   20
  
  # check if any row of x has at least one missing value
  sum(is.na(apply(df.next[,-c(2,3)], 1, sum)))  # 0
  check.col.na = apply(df.next[,-c(2,3)], 2, sum); length(check.col.na); sum(is.na(check.col.na)) # 18 check!   0
  
  # build the model using the training set
  model.next <- train( QNP_obs ~ ., data = train.data.next, method = "glmnet", 
                  trControl = trainControl("cv", number = 10),
                  tuneLength = 10, family="gaussian", 
                  VerboseIter = TRUE, allowParallel = TRUE, trace.it = TRUE
  )
  
# best tuning parameter
  best.lambda.next = model.next$bestTune$lambda # 0.00352
  best.alpha.next  = model.next$bestTune$alpha  # 0.2
  
# coefficient of the final model. You need to specify the best lambda
  keep.next = coef(model.next$finalModel, s = best.lambda.next)
  keep.data.next = as.data.frame(summary(keep.next))
  keep.name.next = row.names(keep.next)
  pred_coef.next <- keep.name.next[keep.data.next[-1,1]]
  keep0.next = as.array(as.matrix(keep.next)); rownames(keep0.next) = NULL
  
# make predictions on the test data
  predictions.next <- model.next %>% predict(test.data.next)
  
# model performance metrics
  keep1.next = data.frame(
    RMSE.next   = RMSE(predictions, test.data.next[, 'QNP_obs']),
    Rsquare.next = R2(predictions, test.data.next[, 'QNP_obs'])
  )	
  keep1.next #       RMSE   Rsquare
  #1 0.4466176 0.8098096
  
  #	  Identify significant predictors
  sum(is.na(keep.next[,1]))	# 0
  sum(abs(keep.next[,1])>0.01)	# 16
  
  keep.select.next = keep.next[ abs(keep.next[,1]) > 0.03,1]  # <<< need ,1 ow. just numerics
  dim(keep.select.next)	# null
  length(keep.select.next)     # 
  list1.next = names(keep.select.next)      # this returns significant ROI names as well as some p.no.
  length(list1.next) # [1] 32
  
  # significant ROIs
  list1.next

# Subset for lm model
    df.std.lm <- data.frame(QNP_obs = df.std.temp[,1],(df.std[,colnames(df.std) %in% list1.next]))

# Lm model with only significant predictors
	model.lm <- lm( QNP_obs ~ ., data = df.std.lm)
  lm(data.frame(scale(model.lm$model)))
	summary(model.lm)
  dwplot(model.lm, dot_args = list(size = 3, pch = 21, fill = "white")) +
    theme_grey(base_size =20) +
    theme(aspect.ratio = 1, legend.position ="none") +
    ylab("Significant predictors") + xlab("Coeffecient values") +
    geom_vline(xintercept = 0, colour = "black", linetype = 2) +
    ggtitle("Base model")
  
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
