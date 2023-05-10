library("tidyverse")
library("table1")
library("tidyr")
library("ggplot2")

getwd()
setwd("~/Desktop/NeuroPathPredict/Input_qnp_data/data wrangling/")
df <- read.csv("qnp_1113_avg_subregions.csv")

colnames(df)[2] = "MFG_per"
colnames(df)[3] = "SMTG_per"
colnames(df)[4] = "IPL_per"

df$HC_per <- rowMeans(df[,c(5,6,7,8,9,10,11,12)], na.rm=TRUE)
df$AM_per <- rowMeans(df[,c(13,14)], na.rm=TRUE)

summary(df)
df_final <- df %>% drop_na()

write.csv(df_final,"qnp_1113_data_ready.csv", row.names = FALSE)
summary(df_final)

ggplot(df_final, aes(x=HC_per)) + geom_histogram()

