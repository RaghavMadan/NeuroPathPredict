#install.packages("tidyverse")  
library(tidyverse)
library("table1")
library("tidyr")
library("ggplot2")

df <- X2022_11_03_Raghav_Madan
df_sub <- subset(df, select= -c(2,3,5,6,7,9,10,11))
df_f <- subset(df_sub, select= c(1,2,3,4))

#Cleaning and extracting HC CA1 data
df_HC_CA1 <- subset(df_sub, select= c(5,6))
colnames(df_HC_CA1)[1] = "HC_CA1_L_per"
colnames(df_HC_CA1)[2] = "HC_CA1_R_per"
df_HC_CA1$HC_CA1_per <- if_else((!is.na(df_HC_CA1$HC_CA1_R_per) & !is.na(df_HC_CA1$HC_CA1_R_per)), (df_HC_CA1$HC_CA1_L_per +df_HC_CA1$HC_CA1_R_per)/2,0)
df_HC_CA1$HC_CA1_per_1 <- ifelse(is.na(df_HC_CA1$HC_CA1_L_per), df_HC_CA1$HC_CA1_R_per, df_HC_CA1$HC_CA1_L_per)
df_HC_CA1$HC_CA1_per_1 <- ifelse(is.na(df_HC_CA1$HC_CA1_R_per), df_HC_CA1$HC_CA1_L_per, df_HC_CA1$HC_CA1_R_per)
df_HC_CA1$HC_CA1_per<- if_else(is.na(df_HC_CA1$HC_CA1_per), df_HC_CA1$HC_CA1_per_1, df_HC_CA1$HC_CA1_per)
df_HC_CA1$HC_CA1_per <- ifelse(df_HC_CA1$HC_CA1_per == 0, df_HC_CA1$HC_CA1_per_1, df_HC_CA1$HC_CA1_per)
summary(df_HC_CA1)

#Adding HC CA1 data to final df
df_f$HC_CA1_per <- df_HC_CA1$HC_CA1_per

#Cleaning and extracting HC CA2 data
df_HC_CA2 <- subset(df_sub, select= c(7,8))
colnames(df_HC_CA2)[1] = "HC_CA2_L_per"
colnames(df_HC_CA2)[2] = "HC_CA2_R_per"
df_HC_CA2$HC_CA2_per <- if_else((!is.na(df_HC_CA2$HC_CA2_R_per) & !is.na(df_HC_CA2$HC_CA2_R_per)), (df_HC_CA2$HC_CA2_L_per +df_HC_CA2$HC_CA2_R_per)/2,0)
df_HC_CA2$HC_CA2_per_1 <- ifelse(is.na(df_HC_CA2$HC_CA2_L_per), df_HC_CA2$HC_CA2_R_per, df_HC_CA2$HC_CA2_L_per)
df_HC_CA2$HC_CA2_per_1 <- ifelse(is.na(df_HC_CA2$HC_CA2_R_per), df_HC_CA2$HC_CA2_L_per, df_HC_CA2$HC_CA2_R_per)
df_HC_CA2$HC_CA2_per<- if_else(is.na(df_HC_CA2$HC_CA2_per), df_HC_CA2$HC_CA2_per_1, df_HC_CA2$HC_CA2_per)
df_HC_CA2$HC_CA2_per <- ifelse(df_HC_CA2$HC_CA2_per == 0, df_HC_CA2$HC_CA2_per_1, df_HC_CA2$HC_CA2_per)
summary(df_HC_CA2)

#Adding HC CA2 data to final df
df_f$HC_CA2_per <- df_HC_CA2$HC_CA2_per

#Cleaning and extracting HC CA3 data
df_HC_CA3 <- subset(df_sub, select= c(9,10))
colnames(df_HC_CA3)[1] = "HC_CA3_L_per"
colnames(df_HC_CA3)[2] = "HC_CA3_R_per"
df_HC_CA3$HC_CA3_per <- if_else((!is.na(df_HC_CA3$HC_CA3_R_per) & !is.na(df_HC_CA3$HC_CA3_R_per)), (df_HC_CA3$HC_CA3_L_per +df_HC_CA3$HC_CA3_R_per)/2,0)
df_HC_CA3$HC_CA3_per_1 <- ifelse(is.na(df_HC_CA3$HC_CA3_L_per), df_HC_CA3$HC_CA3_R_per, df_HC_CA3$HC_CA3_L_per)
df_HC_CA3$HC_CA3_per_1 <- ifelse(is.na(df_HC_CA3$HC_CA3_R_per), df_HC_CA3$HC_CA3_L_per, df_HC_CA3$HC_CA3_R_per)
df_HC_CA3$HC_CA3_per<- if_else(is.na(df_HC_CA3$HC_CA3_per), df_HC_CA3$HC_CA3_per_1, df_HC_CA3$HC_CA3_per)
df_HC_CA3$HC_CA3_per <- ifelse(df_HC_CA3$HC_CA3_per == 0, df_HC_CA3$HC_CA3_per_1, df_HC_CA3$HC_CA3_per)
summary(df_HC_CA3)

#Adding HC CA3 data to final df
df_f$HC_CA3_per <- df_HC_CA3$HC_CA3_per

#Cleaning and extracting HC CA4 data
df_HC_CA4 <- subset(df_sub, select= c(11,12))
colnames(df_HC_CA4)[1] = "HC_CA4_L_per"
colnames(df_HC_CA4)[2] = "HC_CA4_R_per"
df_HC_CA4$HC_CA4_per <- if_else((!is.na(df_HC_CA4$HC_CA4_R_per) & !is.na(df_HC_CA4$HC_CA4_R_per)), (df_HC_CA4$HC_CA4_L_per +df_HC_CA4$HC_CA4_R_per)/2,0)
df_HC_CA4$HC_CA4_per_1 <- ifelse(is.na(df_HC_CA4$HC_CA4_L_per), df_HC_CA4$HC_CA4_R_per, df_HC_CA4$HC_CA4_L_per)
df_HC_CA4$HC_CA4_per_1 <- ifelse(is.na(df_HC_CA4$HC_CA4_R_per), df_HC_CA4$HC_CA4_L_per, df_HC_CA4$HC_CA4_R_per)
df_HC_CA4$HC_CA4_per<- if_else(is.na(df_HC_CA4$HC_CA4_per), df_HC_CA4$HC_CA4_per_1, df_HC_CA4$HC_CA4_per)
df_HC_CA4$HC_CA4_per <- ifelse(df_HC_CA4$HC_CA4_per == 0, df_HC_CA4$HC_CA4_per_1, df_HC_CA4$HC_CA4_per)
summary(df_HC_CA4)

#Adding HC CA4 data to final df
df_f$HC_CA4_per <- df_HC_CA4$HC_CA4_per

#Cleaning and extracting HC DG data
df_HC_DG <- subset(df_sub, select= c(13,14))
colnames(df_HC_DG)[1] = "HC_DG_L_per"
colnames(df_HC_DG)[2] = "HC_DG_R_per"
df_HC_DG$HC_DG_per <- if_else((!is.na(df_HC_DG$HC_DG_R_per) & !is.na(df_HC_DG$HC_DG_R_per)), (df_HC_DG$HC_DG_L_per +df_HC_DG$HC_DG_R_per)/2,0)
df_HC_DG$HC_DG_per_1 <- ifelse(is.na(df_HC_DG$HC_DG_L_per), df_HC_DG$HC_DG_R_per, df_HC_DG$HC_DG_L_per)
df_HC_DG$HC_DG_per_1 <- ifelse(is.na(df_HC_DG$HC_DG_R_per), df_HC_DG$HC_DG_L_per, df_HC_DG$HC_DG_R_per)
df_HC_DG$HC_DG_per<- if_else(is.na(df_HC_DG$HC_DG_per), df_HC_DG$HC_DG_per_1, df_HC_DG$HC_DG_per)
df_HC_DG$HC_DG_per <- ifelse(df_HC_DG$HC_DG_per == 0, df_HC_DG$HC_DG_per_1, df_HC_DG$HC_DG_per)
summary(df_HC_DG)

#Adding HC DG data to final df
df_f$HC_DG_per <- df_HC_DG$HC_DG_per

#Cleaning and extracting HC EC data
df_HC_EC <- subset(df_sub, select= c(15,16))
colnames(df_HC_EC)[1] = "HC_EC_L_per"
colnames(df_HC_EC)[2] = "HC_EC_R_per"
df_HC_EC$HC_EC_per <- if_else((!is.na(df_HC_EC$HC_EC_R_per) & !is.na(df_HC_EC$HC_EC_R_per)), (df_HC_EC$HC_EC_L_per +df_HC_EC$HC_EC_R_per)/2,0)
df_HC_EC$HC_EC_per_1 <- ifelse(is.na(df_HC_EC$HC_EC_L_per), df_HC_EC$HC_EC_R_per, df_HC_EC$HC_EC_L_per)
df_HC_EC$HC_EC_per_1 <- ifelse(is.na(df_HC_EC$HC_EC_R_per), df_HC_EC$HC_EC_L_per, df_HC_EC$HC_EC_R_per)
df_HC_EC$HC_EC_per<- if_else(is.na(df_HC_EC$HC_EC_per), df_HC_EC$HC_EC_per_1, df_HC_EC$HC_EC_per)
df_HC_EC$HC_EC_per <- ifelse(df_HC_EC$HC_EC_per == 0, df_HC_EC$HC_EC_per_1, df_HC_EC$HC_EC_per)
summary(df_HC_EC)

#Adding HC EC data to final df
df_f$HC_EC_per <- df_HC_EC$HC_EC_per

#Cleaning and extracting HC SB data
df_HC_SB <- subset(df_sub, select= c(17,18))
colnames(df_HC_SB)[1] = "HC_SB_L_per"
colnames(df_HC_SB)[2] = "HC_SB_R_per"
df_HC_SB$HC_SB_per <- if_else((!is.na(df_HC_SB$HC_SB_R_per) & !is.na(df_HC_SB$HC_SB_R_per)), (df_HC_SB$HC_SB_L_per +df_HC_SB$HC_SB_R_per)/2,0)
df_HC_SB$HC_SB_per_1 <- ifelse(is.na(df_HC_SB$HC_SB_L_per), df_HC_SB$HC_SB_R_per, df_HC_SB$HC_SB_L_per)
df_HC_SB$HC_SB_per_1 <- ifelse(is.na(df_HC_SB$HC_SB_R_per), df_HC_SB$HC_SB_L_per, df_HC_SB$HC_SB_R_per)
df_HC_SB$HC_SB_per<- if_else(is.na(df_HC_SB$HC_SB_per), df_HC_SB$HC_SB_per_1, df_HC_SB$HC_SB_per)
df_HC_SB$HC_SB_per <- ifelse(df_HC_SB$HC_SB_per == 0, df_HC_SB$HC_SB_per_1, df_HC_SB$HC_SB_per)
summary(df_HC_SB)

#Adding HC SB data to final df
df_f$HC_SB_per <- df_HC_SB$HC_SB_per

#Cleaning and extracting HC TEC data
df_HC_TEC <- subset(df_sub, select= c(19,20))
colnames(df_HC_TEC)[1] = "HC_TEC_L_per"
colnames(df_HC_TEC)[2] = "HC_TEC_R_per"
df_HC_TEC$HC_TEC_per <- if_else((!is.na(df_HC_TEC$HC_TEC_R_per) & !is.na(df_HC_TEC$HC_TEC_R_per)), (df_HC_TEC$HC_TEC_L_per +df_HC_TEC$HC_TEC_R_per)/2,0)
df_HC_TEC$HC_TEC_per_1 <- ifelse(is.na(df_HC_TEC$HC_TEC_L_per), df_HC_TEC$HC_TEC_R_per, df_HC_TEC$HC_TEC_L_per)
df_HC_TEC$HC_TEC_per_1 <- ifelse(is.na(df_HC_TEC$HC_TEC_R_per), df_HC_TEC$HC_TEC_L_per, df_HC_TEC$HC_TEC_R_per)
df_HC_TEC$HC_TEC_per<- if_else(is.na(df_HC_TEC$HC_TEC_per), df_HC_TEC$HC_TEC_per_1, df_HC_TEC$HC_TEC_per)
df_HC_TEC$HC_TEC_per <- ifelse(df_HC_TEC$HC_TEC_per == 0, df_HC_TEC$HC_TEC_per_1, df_HC_TEC$HC_TEC_per)
summary(df_HC_TEC)

#Adding HC TEC data to final df
df_f$HC_TEC_per <- df_HC_TEC$HC_TEC_per

#Cleaning and extracting AM C data
df_AM_C <- subset(df_sub, select= c(21,22))
colnames(df_AM_C)[1] = "AM_C_L_per"
colnames(df_AM_C)[2] = "AM_C_R_per"
df_AM_C$AM_C_per <- if_else((!is.na(df_AM_C$AM_C_R_per) & !is.na(df_AM_C$AM_C_R_per)), (df_AM_C$AM_C_L_per +df_AM_C$AM_C_R_per)/2,0)
df_AM_C$AM_C_per_1 <- ifelse(is.na(df_AM_C$AM_C_L_per), df_AM_C$AM_C_R_per, df_AM_C$AM_C_L_per)
df_AM_C$AM_C_per_1 <- ifelse(is.na(df_AM_C$AM_C_R_per), df_AM_C$AM_C_L_per, df_AM_C$AM_C_R_per)
df_AM_C$AM_C_per<- if_else(is.na(df_AM_C$AM_C_per), df_AM_C$AM_C_per_1, df_AM_C$AM_C_per)
df_AM_C$AM_C_per <- ifelse(df_AM_C$AM_C_per == 0, df_AM_C$AM_C_per_1, df_AM_C$AM_C_per)
summary(df_AM_C)

#Adding AM C data to final df
df_f$AM_C_per <- df_AM_C$AM_C_per

#Cleaning and extracting AM PN data
df_AM_PN <- subset(df_sub, select= c(23,24))
colnames(df_AM_PN)[1] = "AM_PN_L_per"
colnames(df_AM_PN)[2] = "AM_PN_R_per"
df_AM_PN$AM_PN_per <- if_else((!is.na(df_AM_PN$AM_PN_L_per) & !is.na(df_AM_PN$AM_PN_R_per)), (df_AM_PN$AM_PN_L_per +df_AM_PN$AM_PN_R_per)/2,0)
df_AM_PN$AM_PN_per_1 <- ifelse(is.na(df_AM_PN$AM_PN_L_per), df_AM_PN$AM_PN_R_per, df_AM_PN$AM_PN_L_per)
df_AM_PN$AM_PN_per_1 <- ifelse(is.na(df_AM_PN$AM_PN_R_per), df_AM_PN$AM_PN_L_per, df_AM_PN$AM_PN_R_per)
df_AM_PN$AM_PN_per<- if_else(is.na(df_AM_PN$AM_PN_per), df_AM_PN$AM_PN_per_1, df_AM_PN$AM_PN_per)
df_AM_PN$AM_PN_per <- ifelse(df_AM_PN$AM_PN_per == 0, df_AM_PN$AM_PN_per_1, df_AM_PN$AM_PN_per)
summary(df_AM_PN)

#Adding AM PN data to final df
df_f$AM_PN_per <- df_AM_PN$AM_PN_per

write.csv(df_f,"~/Desktop/NeuroPathPredict/NeuroPathPredict/Input_qnp_data/qnp_1113_avg_subregions.csv", row.names = FALSE)
