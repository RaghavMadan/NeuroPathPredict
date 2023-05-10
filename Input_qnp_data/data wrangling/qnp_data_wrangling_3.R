library("tidyverse")
library("table1")
library("tidyr")
library("ggplot2")

getwd()
setwd("~/Desktop/NeuroPathPredict/Input_qnp_data/data wrangling/")
df <- read.csv("qnp_2023102.csv")

df$CA3_L <- rowMeans(df[,c(7,9)], na.rm=TRUE)
df$CA3_R <- rowMeans(df[,c(8,10)], na.rm=TRUE)
df$AM <- rowMeans(df[,c(21,22)], na.rm=TRUE)

colnames(df)[2] = "MFG"
colnames(df)[3] = "SMTG"
colnames(df)[4] = "IPL"
colnames(df)[5] = "CA1_L"
colnames(df)[6] = "CA1_R"
colnames(df)[7] = "CA2_L"
colnames(df)[8] = "CA2_R"
colnames(df)[9] = "CA3_L"
colnames(df)[10] = "CA3_R"
colnames(df)[11] = "CA4_L"
colnames(df)[12] = "CA4_R"
colnames(df)[13] = "DG_L"
colnames(df)[14] = "DG_R"
colnames(df)[15] = "EC_L"
colnames(df)[16] = "EC_R"
colnames(df)[17] = "SB_L"
colnames(df)[18] = "SB_R"
colnames(df)[19] = "TEC_L"
colnames(df)[20] = "TEC_R"
colnames(df)[21] = "AM_L"
colnames(df)[22] = "AM_R"


summary(df)
df_final <- df[c(1,2,3,4,5,6,11,12,13,14,15,16,17,18,23,24,25)]

write.csv(df_final,"qnp_0102_data_ready.csv", row.names = FALSE)
summary(df_final)

ggplot(df_final, aes(x=HC_per)) + geom_histogram()

