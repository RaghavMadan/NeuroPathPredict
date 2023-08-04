## Elastic net regression to identify significant predictors for NPP pipeline
# Code includes all three models (base, r10, r1)
# Y: quantitative tau pathology calculated over each of the 10 regions
# R:4.3.1

library(glmnet) #V:4.1-7
library(Rcpp) #V:1.0.11
library(caret) #V:6.0-94
library(dplyr) #V:1.1.2
library(car)

RNGkind(kind = "Mersenne-Twister")
set.seed(1234)

## Data prep ##
#set directory to your working directory
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
    
# Prepare input data for base, r10, r1 model
    list.cols <- names(df_x)
    
    #Subset for only bz_r1 model
    list.cols.kp.r1 <- list.cols[!grepl("r2|r5|r7|r10|r12.5|r15|r20|r25|r30|r40", list.cols)]
    df.std.bzr1 <- data.frame(QNP_obs = df.std.temp[,1],(df.std[,colnames(df.std) %in% list.cols.kp.r1]))
    dim(df.std.bzr1) # 7610 342
    
    #Subset for only bz_r10 and below model
    list.cols.kp.r10 <- list.cols[!grepl("r12.5|r15|r20|r25|r30|r40", list.cols)]
    df.std.bzr10 <- data.frame(QNP_obs = df.std.temp[,1],(df.std[,colnames(df.std) %in% list.cols.kp.r10]))
    dim(df.std.bzr10) #7610  664

# set up training/testing data for all 3 models
    training.samples.base = defineTrain(nrow(df),0.7)  # get 70% of df as training set # n=7610
    train.data.base = df.std[training.samples.base==0,]; dim(train.data.base) # 5375  931
    test.data.base  = df.std[training.samples.base==1,]; dim(test.data.base) # 2235  931
    
    training.samples.r1 = defineTrain(nrow(df),0.7)
    train.data.r1 = df.std.bzr1[training.samples.r1==0,]; dim(train.data.r1) # 5290  342
    test.data.r1  = df.std.bzr1[training.samples.r1==1,]; dim(test.data.r1) # 2320  342
    
    training.samples.r10 = defineTrain(nrow(df),0.7)
    train.data.r10 = df.std.bzr10[training.samples.r10==0,]; dim(train.data.r10) # 5290  664
    test.data.r10  = df.std.bzr10[training.samples.r10==1,]; dim(test.data.r10) # 2320  664
    
# check if any row of x has at least one missing value
  sum(is.na(apply(df.std[,-c(2,3)], 1, sum)))  # 0
  check.col.na = apply(df.std[,-c(2,3)], 2, sum); length(check.col.na); sum(is.na(check.col.na)) # 929 check!   0

# build the model using the training set
  
  # base model
	model.base <- train( QNP_obs ~ ., data = train.data.base, method = "glmnet", 
           	        trControl = trainControl("cv", number = 10),
                    tuneLength = 10, family="gaussian", 
                    VerboseIter = TRUE, allowParallel = TRUE, trace.it = TRUE
	                )
	
	# r1 model
	model.r1 <- train( QNP_obs ~ ., data = train.data.r1, method = "glmnet", 
	                trControl = trainControl("cv", number = 10),
	                tuneLength = 10, family="gaussian", 
	                VerboseIter = TRUE, allowParallel = TRUE, trace.it = TRUE
	)

	# r10 model
	model.r10 <- train( QNP_obs ~ ., data = train.data.r10, method = "glmnet", 
	                trControl = trainControl("cv", number = 10),
	                tuneLength = 10, family="gaussian", 
	                VerboseIter = TRUE, allowParallel = TRUE, trace.it = TRUE
	)
	
	#par(mfrow=c(3,1))
	#plot(model.base)
	#plot(model.r1)
	#plot(model.r10)
	
# best tuning parameter
	
	param.comp = data.frame(
    lambda.base = model.base$bestTune$lambda,
    lambda.r1 = model.r1$bestTune$lambda, 
    lambda.r10 = model.r10$bestTune$lambda,
    alpha.base  = model.base$bestTune$alpha,
    alpha.r1  = model.r1$bestTune$alpha,
    alpha.r10  = model.r10$bestTune$alpha
	  )
	param.comp
	# lambda.base   lambda.r1  lambda.r10 alpha.base alpha.r1 alpha.r10
	# 0.008278137 0.008292123 0.008119082        0.2      0.2       0.2
  

# make predictions on the test data
	  predictions.base <- model.base %>% predict(test.data.base)
	  predictions.r1 <- model.r1 %>% predict(test.data.r1)
	  predictions.r10 <- model.r10 %>% predict(test.data.r10)
	  
# model performance metrics
	  model.per.metric = data.frame(
		    RMSE.base   = RMSE( predictions.base, test.data.base[, 'QNP_obs']),
		    RMSE.r1     = RMSE( predictions.r1, test.data.r1[, 'QNP_obs']),
		    RMSE.r10    = RMSE( predictions.r10, test.data.r10[, 'QNP_obs']),
		    R2.base     = R2( predictions.base, test.data.base[, 'QNP_obs']),
		    R2.r1       = R2( predictions.r1, test.data.r1[, 'QNP_obs']),
		    R2.r10      = R2( predictions.r10, test.data.r10[, 'QNP_obs'])
		  )	
	  model.per.metric
	  # RMSE.base   RMSE.r1  RMSE.r10   R2.base     R2.r1    R2.r10
	  # 0.4518955 0.4668774 0.4573328 0.8017667 0.7817699 0.7966102

# coefficient of the final model. You need to specify the best lambda
	  lambda.base = model.base$bestTune$lambda
	  lambda.r1 = model.r1$bestTune$lambda
	  lambda.r10 = model.r10$bestTune$lambda

	 #base model
	  keep.base = coef(model.base$finalModel, s = lambda.base)
	  keep.data.base = as.data.frame(summary(keep.base))
	  keep.name.base = row.names(keep.base)
	  pred_coef.base <- keep.name.base[keep.data.base[-1,1]]
	  keep0.base = as.array(as.matrix(keep.base)); rownames(keep0.base) = NULL
	  
	 #r1 model
	  keep.r1 = coef(model.r1$finalModel, s = lambda.r1)
	  keep.data.r1 = as.data.frame(summary(keep.r1))
	  keep.name.r1 = row.names(keep.r1)
	  pred_coef.r1 <- keep.name.r1[keep.data.r1[-1,1]]
	  keep0.r1 = as.array(as.matrix(keep.r1)); rownames(keep0.r1) = NULL
	  
	 #r10 model
	  keep.r10 = coef(model.r10$finalModel, s = lambda.r10)
	  keep.data.r10 = as.data.frame(summary(keep.r10))
	  keep.name.r10 = row.names(keep.r10)
	  pred_coef.r10 <- keep.name.r10[keep.data.r10[-1,1]]
	  keep0.r10 = as.array(as.matrix(keep.r10)); rownames(keep0.r10) = NULL
	  	  
#	  Identify significant predictors
	  
	 #base model
	  sum(is.na(keep.base[,1]))	# 0
	  sum(abs(keep.base[,1])>0.01)	# 747
	  keep.select.base = keep.base[ abs(keep.base[,1]) > 0.01,1]  # <<< thr = 0.01
	  dim(keep.select.base)	# null
	  length(keep.select.base)     # 747
	  list1.base = names(keep.select.base)      # this returns significant ROI names as well as some p.no.
	  list2.base = keep.select.base[!grepl("p.no.", list1.base)]; length(list2.base) # 30
	    # check  
	    sum(grepl("^p.no.", list1.base)) # 717   -  this counts number of entries starting with patient numbers
	    length(list2.base)               #  30   -  this counts non-patient numbers from significant predictors
	    length(keep.select.base)         # 747   -  717 + 30 = 747  checked!

	 #r1 model 
	  sum(is.na(keep.r1[,1]))	# 0
	  sum(abs(keep.r1[,1])>0.01)	# 749
	  keep.select.r1 = keep.r1[ abs(keep.r1[,1]) > 0.01,1]  # <<< thr = 0.01
	  dim(keep.select.r1)	# null
	  length(keep.select.r1)     # 749
	  list1.r1 = names(keep.select.r1)      # this returns significant ROI names as well as some p.no.
	  list2.r1 = keep.select.r1[!grepl("p.no.", list1.r1)]; length(list2.r1) # 33
	    # check  
	      sum(grepl("^p.no.", list1.r1)) # 716   -  this counts number of entries starting with patient numbers
	      length(list2.r1)               #  33   -  this counts non-patient numbers from significant predictors
	      length(keep.select.r1)         # 749   -  716 + 33 = 749  checked!
	      
	 #r10 model 
	  sum(is.na(keep.r10[,1]))	# 0
	  sum(abs(keep.r10[,1])>0.01)	# 733
	  keep.select.r10 = keep.r10[ abs(keep.r10[,1]) > 0.01,1]  # <<< thr = 0.01
	  dim(keep.select.r10)	# null
	  length(keep.select.r10)     # 733
	  list1.r10 = names(keep.select.r10)      # this returns significant ROI names as well as some p.no.
	  list2.r10 = keep.select.r10[!grepl("p.no.", list1.r10)]; length(list2.r10) # 30
	    # check  
	      sum(grepl("^p.no.", list1.r10)) # 703  -  this counts number of entries starting with patient numbers
	      length(list2.r10)               #  30   -  this counts non-patient numbers from significant predictors
	      length(keep.select.r10)         # 733   -  703 + 33 = 733  checked!
	          
 # significant predictors

	 sig.pred <- list(base = names(list2.base), r1 = names(list2.r1), r10 = names(list2.r10))
	 sig.pred
	 # $base
	 # [1] "roiCA4"                   "bz_Buckner2011_3_bz_r10"  "bz_CerebrA_18_bz_r5"     
	 # [4] "bz_CerebrA_23_bz_r10"     "bz_CerebrA_5_bz_r1"       "bz_F_R_bz_r1"            
	 # [7] "bz_MdLF_R_bz_r10"         "bz_Yeo2011_12_bz_r40"     "bz_Yeo2011_6_bz_r40"     
	 # [10] "bz_Yeo2011_9_bz_r7"       "bz_Yeo2011_9_bz_r10"      "bz_Yeo2011_9_bz_r12.5"   
	 # [13] "bz_Yeo2011_9_bz_r15"      "bz_Yeo2011_9_bz_r20"      "bz_Yeo2011_9_bz_r25"     
	 # [16] "bz_Yeo2011_9_bz_r30"      "bz_Yeo2011_9_bz_r40"      "bz_vesRad_thr_0.7_bz_r30"
	 # [19] "bz_vesRad_thr_0.7_bz_r40" "bz_vesRad_thr_1.0_bz_r40" "edt_dxFrom_CerebrA_15"   
	 # [22] "edt_dxFrom_CerebrA_25"    "edt_dxFrom_CerebrA_3"     "edt_dxFrom_CerebrA_60"   
	 # [25] "edt_dxFrom_CerebrA_7"     "edt_dxFrom_CerebrA_82"    "edt_dxFrom_Yeo2011_10"   
	 # [28] "edt_dxFrom_Yeo2011_3"     "var_CerebrA_5"            "var_F_R"                 
	 # 
	 # $r1
	 # [1] "(Intercept)"           "roiCA4"                "roiDG"                 "bz_CST_R_bz_r1"       
	 # [5] "bz_C_PH_R_bz_r1"       "bz_CerebrA_18_bz_r1"   "bz_Yeo2011_7_bz_r1"    "edt_dxFrom_CNIII_R"   
	 # [9] "edt_dxFrom_C_FP_L"     "edt_dxFrom_CerebrA_15" "edt_dxFrom_CerebrA_25" "edt_dxFrom_CerebrA_36"
	 # [13] "edt_dxFrom_CerebrA_3"  "edt_dxFrom_CerebrA_45" "edt_dxFrom_CerebrA_60" "edt_dxFrom_CerebrA_67"
	 # [17] "edt_dxFrom_CerebrA_7"  "edt_dxFrom_CerebrA_82" "edt_dxFrom_CerebrA_84" "edt_dxFrom_CerebrA_98"
	 # [21] "edt_dxFrom_ILF_R"      "edt_dxFrom_SLF1_L"     "edt_dxFrom_SLF2_L"     "edt_dxFrom_TR_S_L"    
	 # [25] "edt_dxFrom_Yeo2011_10" "edt_dxFrom_Yeo2011_11" "edt_dxFrom_Yeo2011_3"  "edt_dxFrom_Yeo2011_6" 
	 # [29] "var_CPT_O_R"           "var_CPT_P_R"           "var_CST_R"             "var_CerebrA_5"        
	 # [33] "var_Yeo2011_7"        
	 # 
	 # $r10
	 # [1] "roiCA4"                "bz_C_PH_R_bz_r2"       "bz_CerebrA_18_bz_r5"   "bz_CerebrA_18_bz_r7"  
	 # [5] "bz_CerebrA_23_bz_r10"  "bz_CerebrA_36_bz_r10"  "bz_CerebrA_5_bz_r1"    "bz_MdLF_R_bz_r10"     
	 # [9] "bz_Yeo2011_9_bz_r7"    "bz_Yeo2011_9_bz_r10"   "edt_dxFrom_CNIII_L"    "edt_dxFrom_CNIII_R"   
	 # [13] "edt_dxFrom_CerebrA_15" "edt_dxFrom_CerebrA_25" "edt_dxFrom_CerebrA_36" "edt_dxFrom_CerebrA_3" 
	 # [17] "edt_dxFrom_CerebrA_4"  "edt_dxFrom_CerebrA_60" "edt_dxFrom_CerebrA_67" "edt_dxFrom_CerebrA_7" 
	 # [21] "edt_dxFrom_CerebrA_82" "edt_dxFrom_CerebrA_84" "edt_dxFrom_MCP"        "edt_dxFrom_TR_S_L"    
	 # [25] "edt_dxFrom_Yeo2011_10" "edt_dxFrom_Yeo2011_3"  "edt_dxFrom_Yeo2011_6"  "edt_dxFrom_Yeo2011_9" 
	 # [29] "var_CerebrA_5"         "var_F_R" 

#### Exploratory factor analysis ######
	 ##Manual removal of bz predictors with larger radius if more than one is picked
	 #base model
	  sig.cov.rm.base <- c("roiCA4","bz_Yeo2011_9_bz_r10","bz_Yeo2011_9_bz_r12.5",
	                  "bz_Yeo2011_9_bz_r15", "bz_Yeo2011_9_bz_r20", "bz_Yeo2011_9_bz_r25",
	                  "bz_Yeo2011_9_bz_r30","bz_Yeo2011_9_bz_r40","bz_vesRad_thr_0.7_bz_r40")
	  sig.cov.base = names(list2.base)[!(names(list2.base) %in% sig.cov.rm.base)]
	  length(sig.cov.base) #21
	  sig.cov.base
	  
	  #r10 model
	  sig.cov.rm.r10 <- c("roiCA4","bz_CerebrA_18_bz_r7", "bz_Yeo2011_9_bz_r10")
	  sig.cov.r10 = names(list2.r10)[!(names(list2.r10) %in% sig.cov.rm.r10)]
	  length(sig.cov.r10) #27
	  sig.cov.r10
	  
	 # data prep
	  # base model
	  df_x.next.base <- df.std[,colnames(df.std) %in% sig.cov.base]
	  df.next.prep.base <- data.frame(QNP_obs = df.std$QNP_obs, roi = df.std$roi, df_x.next.base)
	  
	  df.next.base <- data.frame(model.matrix( ~ ., df.next.prep.base, ignore.intercept = TRUE)) %>%
	    select(-c(roiMFG, roiSMTG, roiIPL, roiCA1, roiCA3, roiEC, roiSB, roiDG,X.Intercept.))
	  
	  # r1 model
	  df_x.next.r1 <- df.std[,colnames(df.std) %in% names(list2.r1)]
	  df.next.prep.r1 <- data.frame(QNP_obs = df.std$QNP_obs, roi = df.std$roi, df_x.next.r1)
	  
	  df.next.r1 <- data.frame(model.matrix( ~ ., df.next.prep.r1, ignore.intercept = TRUE)) %>%
	    select(-c(roiMFG, roiSMTG, roiIPL, roiCA1, roiCA3, roiEC, roiSB, roiDG,X.Intercept.))
	  
	  # Lm model with only significant predictors
	    #base model
	    model.lm.base <- lm( QNP_obs ~ ., data = df.next.base)
	    lm(data.frame(scale(model.lm.base$model)))
	    summary(model.lm.base)
	    
	     # linearly dependent variables
	      ld.vars.base <- attributes(alias(model.lm.base)$Complete)$dimnames[[1]]
	      ld.vars.base
	      
	      df.next.base.ld.vars.prep <- df.next.base[,colnames(df.next.base) %in% ld.vars.base]
	      df.next.base.ld.vars <- data.frame(QNP_obs = df.std$QNP_obs, df.next.base.ld.vars.prep)
	  
	    model.lm.base.ld.vars <- lm( QNP_obs ~., data=df.next.base.ld.vars)
	    summary(model.lm.base.ld.vars)
	    
	   #r1 model
	    model.lm.r1 <- lm( QNP_obs ~ ., data = df.next.r1)
	    lm(data.frame(scale(model.lm.r1$model)))
	    summary(model.lm.r1)
	    
	    # linearly dependent variables r1 model
	    ld.vars.r1 <- attributes(alias(model.lm.r1)$Complete)$dimnames[[1]]
	    ld.vars.r1
	    
	    df.next.r1.ld.vars.prep <- df.next.r1[,colnames(df.next.r1) %in% ld.vars.r1]
	    df.next.r1.ld.vars <- data.frame(QNP_obs = df.std$QNP_obs, df.next.r1.ld.vars.prep)
	    
	    model.lm.r1.ld.vars <- lm( QNP_obs ~., data=df.next.r1.ld.vars)
	    summary(model.lm.r1.ld.vars)

	    ###### stopped here since the results of all lm model have adjusted R2 = 0.455 #####
	    ## Did not run lm for R10 model ##
	  

	  