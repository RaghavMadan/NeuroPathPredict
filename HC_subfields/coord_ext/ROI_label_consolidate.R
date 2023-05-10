

library("table1")
library(tidyverse)

df <- read.table("~/Desktop/NeuroPathPredict/HC_subfields/coord_ext/ROI_HC_xyz_coord_labels.csv",sep = ",", header = TRUE)
df <- df[-c(5)]
vars <-c("X","X0","X1","X2")
df_hc <- df[vars]

df_hc$roi_label <- paste(df$roi_label_y, df$roi_label_x.1, df$roi_label_y.1, df$roi_label_x.2, df$roi_label_y.2,
                         df$roi_label_x.3, df$roi_label_y.3, df$roi_label_x.4, df$roi_label_y.4, 
                         df$roi_label_x.5, df$roi_label_y.5, df$roi_label_x.6, df$roi_label_y.6, df$roi_label)

df_hc <- data.frame(lapply(df_hc, trimws), stringsAsFactors = FALSE)
df_hc$X0 <- as.numeric(df_hc$X0)
df_hc$X1 <- as.numeric(df_hc$X1)
df_hc$X2 <- as.numeric(df_hc$X2)
df_hc$X <- as.integer(df_hc$X)

write_csv(df_hc,"~/Desktop/NeuroPathPredict/HC_subfields/coord_ext/ROI_HC_subfield_labels.csv")
