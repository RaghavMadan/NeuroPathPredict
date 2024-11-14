## Elastic net regression to identify significant predictors for NPP pipeline V1 model
# Y: quantitative tau pathology calculated for 10 participants for region MFG
# R:4.3.1
# Base model
# Aliased variables have also been removed compared to last version. 


# ???? Check versions ????
library(glmnet) #V:4.1-7
library(Rcpp) #V:1.0.10
library(caret) #V:6.0-94
library(dplyr) #V:1.1.2
library(dotwhisker) #V:0.7.4
library(data.table)
library(ggplot2)
library(doParallel)
library(GGally)
library(ggcorrplot)
library(broom)
library(reshape2)

## Data prep ##

RNGkind(kind = "Mersenne-Twister")
set.seed(1234)

setwd("~/Desktop/NeuroPathPredict/V1/Model_V1/Input_data")

# define standardize()

#standardize = function(x){
#  tm = mean(x, na.rm=T)
#  ts = sd(x, na.rm=T)
#  temp = (x-tm)/ts
#  return(temp)
#}

#df.std.temp = apply(df_in[,-c(1:3)], 2, function(x) standardize(x))
#df.std <- data.frame(QNP_obs = df.std.temp[,1], sub.id = df_in[,3], df.std.temp[,-c(1)])
#dim(df.std)  # 438571    450

### Load standardized and processed df ###
df.std <- read.csv("XY_MFG_V1_right_std_1021.csv")

# Prepare group variables for predictors
df.std$sub.id <- as.factor(df.std$sub.id)  # Convert sub.id to factor
class(df.std$sub.id)  # "factor"

# Data check
names(df.std)[1:20]
dim(df.std)	#438571    450

### Data visualization before implementing model ###

# Create a data frame for plotting
plot_data <- data.frame(sub.id = df.std$sub.id, QNP_obs = df.std$QNP_obs)

# Plot QNP_obs by sub.id using box plots
ggplot(plot_data, aes(x = sub.id, y = QNP_obs)) +
  geom_boxplot(fill = "blue", alpha = 0.7) +
  labs(title = "QNP_obs by sub.id",
       x = "sub.id",
       y = "QNP_obs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Calculate quartiles and deciles
#quartiles <- quantile(QNP_obs, probs = seq(0, 1, 0.25))
#deciles <- quantile(QNP_obs, probs = seq(0, 1, 0.1))

# Categorize QNP_obs into quartiles and deciles
#plot_data$Quartile <- cut(QNP_obs, breaks = quartiles, include.lowest = TRUE, labels = FALSE)
#plot_data$Decile <- cut(QNP_obs, breaks = deciles, include.lowest = TRUE, labels = FALSE)

# Plot the distribution within each quartile
#ggplot(plot_data, aes(x = factor(Quartile), y = QNP_obs)) +
#  geom_boxplot(fill = "blue", alpha = 0.7) +
#  labs(title = "Distribution of QNP_obs within Quartiles",
#       x = "Quartiles",
#       y = "QNP_obs") +
#  theme_minimal()

# Plot the distribution within each decile
#ggplot(plot_data, aes(x = factor(Decile), y = QNP_obs)) +
#  geom_boxplot(fill = "green", alpha = 0.7) +
#  labs(title = "Distribution of QNP_obs within Deciles",
#       x = "Deciles",
#       y = "QNP_obs") +
#  theme_minimal()

#df.std <- as.data.table(df.std)

# define function to split data into training and test sets
# n: number of observations
# p: proportion of observations in training set

defineTrain = function(n,p){
  choice1 = runif(n,0,1)
  choice2 = rep(0, n)
  choice2[choice1>p] = 1
  return(choice2)
}

# set up training/testing data
training.samples = defineTrain(nrow(df.std),0.8)  # get 80% of df as training set # n=438571
train.data = df.std[training.samples==0,]; dim(train.data) # 351070    450
test.data  = df.std[training.samples==1,]; dim(test.data) # 87501   450

# check if any row of x has at least one missing value
sum(is.na(apply(df.std[,-c(1,2,3,4)], 1, sum)))  # 0
check.col.na = apply(df.std[,-c(1,2,3)], 2, sum); length(check.col.na); sum(is.na(check.col.na)) # 447 check!   0

# Call garbage collection
gc()

# Build the model using the training set
model <- train(QNP_obs ~ ., data = train.data, method = "glmnet", 
               trControl = trainControl("cv", number = 10),
               tuneLength = 10, family = "gaussian", 
               verboseIter = TRUE, allowParallel = TRUE, trace.it = TRUE)

# best tuning parameter
best.lambda = model$bestTune$lambda # # 0.000146
best.alpha  = model$bestTune$alpha  # 0.6

# coefficient of the final model. You need to specify the best lambda
keep = coef(model$finalModel, s = best.lambda)
keep.data = as.data.frame(summary(keep))
keep.name = row.names(keep)
pred_coef <- keep.name[keep.data[-1,1]]
keep0 = as.array(as.matrix(keep)); rownames(keep0) = NULL

# make predictions on the test data
predictions <- model %>% predict(test.data)

# Convert to numeric if necessary
test.data[, 'QNP_obs'] <- as.numeric(unlist(test.data[, 'QNP_obs']))

# Convert test.data[, 'QNP_obs'] to a plain numeric vector
QNP_obs_numeric <- as.numeric(unlist(test.data[, 'QNP_obs']))

# model performance metrics
keep1 = data.frame(
  RMSE    = RMSE(predictions, QNP_obs_numeric),
  Rsquare = R2(predictions, QNP_obs_numeric)
)
keep1
#      RMSE   Rsquare
# 0.760606 0.4341046

#	  Identify significant predictors
sum(is.na(keep[,1]))	# 0
sum(abs(keep[,1])>0.01)	# 272

keep.select = keep[ abs(keep[,1]) > 0.9,1]  # <<< thr = 0.9
dim(keep.select)	# null
length(keep.select)     # 30
list1 = names(keep.select)      # this returns significant ROI names as well as some p.no.
list2 = keep.select[!grepl("sub.id", list1)]; length(list2) # 27

# check  
sum(grepl("^sub.id.", list1)) # 0   -  this counts number of entries starting with patient numbers
length(list2)               #  30   -  this counts non-patient numbers from significant predictors
length(keep.select)         # 30   

# significant ROIs
names(list2)

# Save the names of list2 to a text file
write(names(list2), file = "EN1_SigPred_30.txt")

#### Exploratory factor analysis ######
sig.cov1 = names(list2)
length(sig.cov1) #30
sig.cov1

df_x.next <- df.std[, colnames(df.std) %in% sig.cov1]
df_next <- data.frame(QNP_obs = df.std$QNP_obs, df_x.next)

# Lm model 1 with only 30 significant predictors
model.lm <- lm( QNP_obs ~ ., data = df_next)
summary(model.lm)

  
##################### 2nd Elastic Net with significant predictors only ###################

# set up training/testing data
training.samples_next = defineTrain(nrow(df),0.8)  # get 80% of df as training set # n=438571
train.data_next = df_next[training.samples_next==0,]; dim(train.data_next) # 350389    31
test.data_next  = df_next[training.samples_next==1,]; dim(test.data_next) # 88182   31

# Call garbage collection
gc()

# Build the model using the training set
model_next <- train(QNP_obs ~ ., data = train.data_next, method = "glmnet", 
               trControl = trainControl("cv", number = 10),
               tuneLength = 10, family = "gaussian", 
               verboseIter = TRUE, allowParallel = TRUE, trace.it = TRUE)

# best tuning parameter
best.lambda_next = model_next$bestTune$lambda # # 0.0001374
best.alpha_next  = model_next$bestTune$alpha  # 1

# coefficient of the final model. You need to specify the best lambda
keep_next = coef(model_next$finalModel, s = best.lambda_next)
keep.data_next = as.data.frame(summary(keep_next))
keep.name_next = row.names(keep_next)
pred_coef_next <- keep.name_next[keep.data_next[-1,1]]
keep0_next = as.array(as.matrix(keep_next)); rownames(keep0_next) = NULL

# make predictions on the test data
predictions_next <- model_next %>% predict(test.data_next)

# Convert to numeric if necessary
test.data_next[, 'QNP_obs'] <- as.numeric(unlist(test.data_next[, 'QNP_obs']))

# Convert test.data[, 'QNP_obs'] to a plain numeric vector
QNP_obs_numeric_next <- as.numeric(unlist(test.data_next[, 'QNP_obs']))

# model performance metrics
keep1_next = data.frame(
  RMSE    = RMSE(predictions_next, QNP_obs_numeric_next),
  Rsquare = R2(predictions_next, QNP_obs_numeric_next)
)
keep1_next
#      RMSE   Rsquare
# 0.8872191 0.2218573

#	  Identify significant predictors
sum(is.na(keep_next[,1]))	# 0
sum(abs(keep_next[,1])>0.01)	# 29

keep.select_next = keep_next[ abs(keep_next[,1]) > 0.3,1]  # <<< thr = 0.3
dim(keep.select_next)	# null
length(keep.select_next)     # 20
list1_next = names(keep.select_next)      # this returns significant ROI names as well as some p.no.
list2_next = keep.select_next[!grepl("sub.id", list1_next)]; length(list2_next) # 14

# check  
length(list2_next)               #  20   -  this counts non-patient numbers from significant predictors
length(keep.select_next)         # 20   

# significant ROIs
names(list2_next)

sig.cov1_next = names(list2_next)

df_x.next2 <- df_next[, colnames(df_next) %in% sig.cov1_next]
df_next2 <- data.frame(QNP_obs = df_next$QNP_obs, df_x.next2)

# Lm model 2 with 20 significant predictors
model.lm2 <- lm( QNP_obs ~ ., data = df_next2)
summary(model.lm2)
