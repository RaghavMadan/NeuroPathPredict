#install.packages("glmnet")
#install.packages('Rcpp')
#install.packages("MASS")

library('Rcpp')
library("glmnet")
library("mice")
library("Hmisc")
library("GeneNet")
library("MASS")
library("car")

getwd()
setwd("~/Desktop/NeuroPathPredict/Model_V0/")
df_y<- read.csv("Y_qnp_data_0426.csv")
df_x<- read.csv("X_cov_roi_0515.csv")

rownames(df_x) <- df_x$rois
rownames(df_y) <- df_y$UWA.

df_y<- df_y[-c(1)]
df_x<- df_x[-c(1)]

######## Process X ########
# z-score standardization for columns 3:931
df_x_scaled <- df_x
df_x_scaled[,3:931] <- scale(df_x[,3:931])

# check standard variation and mean of each column and store in data frame
X <- data.frame(colnames(df_x_scaled[,3:931]))
X$mean <- apply(df_x_scaled[,3:931], 2, mean)
X$sd <- apply(df_x_scaled[,3:931], 2, sd)
X$var <- apply(df_x_scaled[,3:931], 2, var)


write.csv(df_x_scaled, "X_cov_roi_0524.csv", row.names = TRUE)


#x_MFG <- cbind(df_x[1,])
#x_MFG <- x_MFG[rep(1,761),]
#x_SMTG <- cbind(df_x[2,])
#x_SMTG <- x_SMTG[rep(1,761),]
#x_IPL <- cbind(df_x[3,])
#x_IPL <- x_IPL[rep(1,761),]
# x_AM <- cbind(df_x[4,])
# x_AM <- x_AM[rep(1,761),]
# x_CA1 <- cbind(df_x[5,])
# x_CA1 <- x_CA1[rep(1,761),]
# x_CA3 <- cbind(df_x[6,])
# x_CA3 <- x_CA3[rep(1,761),]
# x_CA4 <- cbind(df_x[7,])
# x_CA4 <- x_CA4[rep(1,761),]
# x_DG <- cbind(df_x[8,])
# x_DG <- x_DG[rep(1,761),]
# x_EC <- cbind(df_x[9,])
# x_EC <- x_EC[rep(1,761),]
# x_SB <- cbind(df_x[10,])
# x_SB <- x_SB[rep(1,761),]
# X <- rbind(x_MFG,x_SMTG,x_IPL,x_AM,x_CA1,x_CA3,x_CA4,x_DG,x_EC,x_SB)



####### Process Y #######
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
df_y["CA1",][is.na(df_y["CA1",])] <- a["CA1",]
df_y["CA3",][is.na(df_y["CA3",])] <- a["CA3",]
df_y["CA4",][is.na(df_y["CA4",])] <- a["CA4",]
df_y["DG",][is.na(df_y["DG",])] <- a["DG",]
df_y["SB",][is.na(df_y["SB",])] <- a["SB",]
df_y["EC",][is.na(df_y["EC",])] <- a["EC",]

summary(Y)
# calculate mean, sd, var for Y ##
Y_sum <- data.frame(colnames(Y))
Y_sum$mean <- apply(Y, 2, mean)
Y_sum$sd <- apply(Y, 2, sd)
Y_sum$var <- apply(Y, 2, var)


t_y <- data.frame(t(df_y))
t_ybc <- t_y

##Calculating lambda values for box-cox transformation

b_MFG <- boxcox(lm(t_y$MFG~1))
l_MFG <- b_MFG$x[which.max(b_MFG$y)]
t_ybc$MFG <- (t_y$MFG^ l_MFG -1)/l_MFG

b_SMTG <- boxcox(lm(t_y$SMTG~1))
l_SMTG <- b_SMTG$x[which.max(b_SMTG$y)]
t_ybc$SMTG <- (t_y$SMTG ^ l_SMTG -1)/l_SMTG

b_IPL <- boxcox(lm(t_y$IPL~1))
l_IPL <- b_IPL$x[which.max(b_IPL$y)]
t_ybc$IPL <- (t_y$IPL ^ l_IPL -1)/l_IPL

b_AM <- boxcox(lm(t_y$AM~1))
l_AM <- b_AM$x[which.max(b_AM$y)]
t_ybc$AM <- (t_y$AM ^ l_AM -1)/l_AM

b_CA1 <- boxcox(lm(t_y$CA1~1))
l_CA1 <- b_CA1$x[which.max(b_CA1$y)]
t_ybc$CA1 <- (t_y$CA1 ^ l_CA1 -1)/l_CA1

b_CA3 <- boxcox(lm(t_y$CA3~1))
l_CA3 <- b_CA3$x[which.max(b_CA3$y)]
t_ybc$CA3 <- (t_y$CA3 ^ l_CA3 -1)/l_CA3

b_CA4 <- boxcox(lm(t_y$CA4~1))
l_CA4 <- b_CA4$x[which.max(b_CA4$y)]
t_ybc$CA4 <- (t_y$CA4 ^ l_CA4 -1)/l_CA4

b_DG <- boxcox(lm(t_y$DG~1))
l_DG <- b_DG$x[which.max(b_DG$y)]
t_ybc$DG <- (t_y$DG ^ l_DG -1)/l_DG

b_SB <- boxcox(lm(t_y$SB~1))
l_SB <- b_SB$x[which.max(b_SB$y)]
t_ybc$SB <- (t_y$SB ^ l_SB -1)/l_SB

b_EC <- boxcox(lm(t_y$EC~1))
l_EC <- b_EC$x[which.max(b_EC$y)]
t_ybc$EC <- (t_y$EC ^ l_EC -1)/l_EC

df_ybc <- data.frame(t(t_ybc))

y <- cbind(df_ybc[1,],df_ybc[2,],df_ybc[3,],df_ybc[4,],
df_ybc[5,],df_ybc[6,],df_ybc[7,],df_ybc[8,],df_ybc[9,],df_ybc[10,])
Y <- as.data.frame(t(y))

write.csv(Y, "Y_qnp_data_0524.csv")

##Plotting distibutions (raw, log transform, BC transform)
#MFG

attach(t_y)
par(mfrow=c(1,3))
hist(t_y$MFG, breaks = 50)
qqnorm(t_y$MFG)
qqline(t_y$MFG)
boxplot(t_y$MFG, horizontal = TRUE, main = "Boxplot of MFG")

attach(t_ylog)
par(mfrow=c(1,3))
hist(t_ylog$MFG, breaks = 50)
qqnorm(t_ylog$MFG)
qqline(t_ylog$MFG)
boxplot(t_ylog$MFG, horizontal = TRUE, main = "Boxplot of log(MFG)")

attach(t_ybc)
par(mfrow=c(1,3))
hist(t_ybc$MFG, breaks = 50)
qqnorm(t_ybc$MFG)
qqline(t_ybc$MFG)
boxplot(t_ybc$MFG, horizontal = TRUE, main = "Boxplot of bc(MFG)")

#SMTG

attach(t_y)
par(mfrow=c(1,3))
hist(t_y$SMTG, breaks = 50)
qqnorm(t_y$SMTG)
qqline(t_y$SMTG)
boxplot(t_y$SMTG, horizontal = TRUE, main = "Boxplot of SMTG")

attach(t_ylog)
par(mfrow=c(1,3))
hist(t_ylog$SMTG, breaks = 50)
qqnorm(t_ylog$SMTG)
qqline(t_ylog$SMTG)
boxplot(t_ylog$SMTG, horizontal = TRUE, main = "Boxplot of log(SMTG)")

attach(t_ybc)
par(mfrow=c(1,3))
hist(t_ybc$SMTG, breaks = 50)
qqnorm(t_ybc$SMTG)
qqline(t_ybc$SMTG)
boxplot(t_ybc$SMTG, horizontal = TRUE, main = "Boxplot of bc(SMTG)")

##IPL

attach(t_y)
par(mfrow=c(1,3))
hist(t_y$IPL, breaks = 50)
qqnorm(t_y$IPL)
qqline(t_y$IPL)
boxplot(t_y$IPL, horizontal = TRUE, main = "Boxplot of IPL")

attach(t_ylog)
par(mfrow=c(1,3))
hist(t_ylog$IPL, breaks = 50)
qqnorm(t_ylog$IPL)
qqline(t_ylog$IPL)
boxplot(t_ylog$IPL, horizontal = TRUE, main = "Boxplot of log(IPL)")

attach(t_ybc)
par(mfrow=c(1,3))
hist(t_ybc$IPL, breaks = 50)
qqnorm(t_ybc$IPL)
qqline(t_ybc$IPL)
boxplot(t_ybc$IPL, horizontal = TRUE, main = "Boxplot of bc(IPL)")

##AM

par(mfrow=c(1,3))
hist(t_y$AM, breaks = 50)
qqnorm(t_y$AM)
qqline(t_y$AM)
boxplot(t_y$AM, horizontal = TRUE, main = "Boxplot of AM")

par(mfrow=c(1,3))
hist(t_ylog$AM, breaks = 50)
qqnorm(t_ylog$AM)
qqline(t_ylog$AM)
boxplot(t_ylog$AM, horizontal = TRUE, main = "Boxplot of log(AM)")

attach(t_ybc)
par(mfrow=c(1,3))
hist(t_ybc$AM, breaks = 50)
qqnorm(t_ybc$AM)
qqline(t_ybc$AM)
boxplot(t_ybc$AM, horizontal = TRUE, main = "Boxplot of bc(AM)")

##CA1

attach(t_y)
par(mfrow=c(1,3))
hist(t_y$CA1, breaks = 50)
qqnorm(t_y$CA1)
qqline(t_y$CA1)
boxplot(t_y$CA1, horizontal = TRUE, main = "Boxplot of CA1")

attach(t_ylog)
par(mfrow=c(1,3))
hist(t_ylog$CA1, breaks = 50)
qqnorm(t_ylog$CA1)
qqline(t_ylog$CA1)
boxplot(t_ylog$CA1, horizontal = TRUE, main = "Boxplot of log(CA1)")

attach(t_ybc)
par(mfrow=c(3,1))
hist(t_ybc$CA1, breaks = 50)
qqnorm(t_ybc$CA1)
qqline(t_ybc$CA1)
boxplot(t_ybc$CA1, horizontal = TRUE, main = "Boxplot of bc(CA1)")


##CA3

attach(t_y)
par(mfrow=c(1,3))
hist(t_y$CA3, breaks = 50)
qqnorm(t_y$CA3)
qqline(t_y$CA3)
boxplot(t_y$CA3, horizontal = TRUE, main = "Boxplot of CA3")

attach(t_ylog)
par(mfrow=c(1,3))
hist(t_ylog$CA3, breaks = 50)
qqnorm(t_ylog$CA3)
qqline(t_ylog$CA3)
boxplot(t_ylog$CA3, horizontal = TRUE, main = "Boxplot of log(CA3)")

attach(t_ybc)
par(mfrow=c(1,3))
hist(t_ybc$CA3, breaks = 50)
qqnorm(t_ybc$CA3)
qqline(t_ybc$CA3)
boxplot(t_ybc$CA3, horizontal = TRUE, main = "Boxplot of bc(CA3)")

##CA4

attach(t_y)
par(mfrow=c(3,1))
hist(t_y$CA4, breaks = 50)
qqnorm(t_y$CA4)
qqline(t_y$CA4)
boxplot(t_y$CA4, horizontal = TRUE, main = "Boxplot of CA4")

attach(t_ylog)
par(mfrow=c(3,1))
hist(t_ylog$CA4, breaks = 50)
qqnorm(t_ylog$CA4)
qqline(t_ylog$CA4)
boxplot(t_ylog$CA4, horizontal = TRUE, main = "Boxplot of log(CA4)")

attach(t_ybc)
par(mfrow=c(3,1))
hist(t_ybc$CA4, breaks = 50)
qqnorm(t_ybc$CA4)
qqline(t_ybc$CA4)
boxplot(t_ybc$CA4, horizontal = TRUE, main = "Boxplot of bc(CA4)")

##DG

attach(t_y)
par(mfrow=c(3,1))
hist(t_y$DG, breaks = 50)
qqnorm(t_y$DG)
qqline(t_y$DG)
boxplot(t_y$DG, horizontal = TRUE, main = "Boxplot of DG")

attach(t_ylog)
par(mfrow=c(3,1))
hist(t_ylog$DG, breaks = 50)
qqnorm(t_ylog$DG)
qqline(t_ylog$DG)
boxplot(t_ylog$DG, horizontal = TRUE, main = "Boxplot of log(DG)")

attach(t_ybc)
par(mfrow=c(3,1))
hist(t_ybc$DG, breaks = 50)
qqnorm(t_ybc$DG)
qqline(t_ybc$DG)
boxplot(t_ybc$DG, horizontal = TRUE, main = "Boxplot of bc(DG)")

##EC

attach(t_y)
par(mfrow=c(3,1))
hist(t_y$EC, breaks = 50)
qqnorm(t_y$EC)
qqline(t_y$EC)
boxplot(t_y$EC, horizontal = TRUE, main = "Boxplot of EC")

attach(t_ylog)
par(mfrow=c(3,1))
hist(t_ylog$EC, breaks = 50)
qqnorm(t_ylog$EC)
qqline(t_ylog$EC)
boxplot(t_ylog$EC, horizontal = TRUE, main = "Boxplot of log(EC)")

attach(t_ybc)
par(mfrow=c(3,1))
hist(t_ybc$EC, breaks = 50)
qqnorm(t_ybc$EC)
qqline(t_ybc$EC)
boxplot(t_ybc$EC, horizontal = TRUE, main = "Boxplot of bc(EC)")

##SB

attach(t_y)
par(mfrow=c(3,1))
hist(t_y$SB, breaks = 50)
qqnorm(t_y$SB)
qqline(t_y$SB)
boxplot(t_y$SB, horizontal = TRUE, main = "Boxplot of SB")

attach(t_ylog)
par(mfrow=c(3,1))
hist(t_ylog$SB, breaks = 50)
qqnorm(t_ylog$SB)
qqline(t_ylog$SB)
boxplot(t_ylog$SB, horizontal = TRUE, main = "Boxplot of log(SB)")

attach(t_ybc)
par(mfrow=c(3,1))
hist(t_ybc$SB, breaks = 50)
qqnorm(t_ybc$SB)
qqline(t_ybc$SB)
boxplot(t_ybc$SB, horizontal = TRUE, main = "Boxplot of bc(SB)")


