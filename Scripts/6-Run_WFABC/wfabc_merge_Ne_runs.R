library(boa)
library (dplyr)

s_H1N1 <- read.table (snakemake@input[[1]])
mean_s <- rowMeans(s_H1N1)
mean_95 <- apply(s_H1N1,1,boa.hpd,1-0.95)
mean_s <- rbind( mean_s, mean_95)
mean_s <- t(mean_s)
mean_s <- mean_s %>% as.data.frame () %>% rename ( mean_lower_95= "Lower Bound", mean_upper_95= "Upper Bound")


minus_sd <- read.table (snakemake@input[[2]])
minus_sd_s <- rowMeans(minus_sd)
minus_95 <- apply(minus_sd,1,boa.hpd,1-0.95)
minus_sd_s <- rbind( minus_sd_s, minus_95)
minus_sd_s <- t(minus_sd_s)
minus_sd_s <- minus_sd_s %>% as.data.frame () %>% rename ( minus_sd_lower_95= "Lower Bound", minus_sd_upper_95= "Upper Bound")

plus_sd <- read.table (snakemake@input[[3]])
plus_sd_s <- rowMeans(plus_sd)
plus_95 <- apply(plus_sd,1,boa.hpd,1-0.95)
plus_sd_s <- rbind( plus_sd_s, minus_95)
plus_sd_s <- t(plus_sd_s)
plus_sd_s <- plus_sd_s %>% as.data.frame () %>% rename ( plus_sd_lower_95= "Lower Bound", plus_sd_upper_95= "Upper Bound")

Names <- read.csv (snakemake@input[[4]])
Names <- select (Names, hhsubid_pos)
selection <- cbind (Names,mean_s,minus_sd_s,plus_sd_s)

selection <- selection %>% mutate ( mean_pos_selection = ifelse (mean_lower_95 >0 &  mean_upper_95>0, "TRUE", "FALSE")) %>% mutate ( minus_sd_pos_selection = ifelse (minus_sd_lower_95 >0 &  minus_sd_upper_95>0, "TRUE", "FALSE"))%>% mutate ( plus_sd_pos_selection = ifelse (plus_sd_lower_95 >0 &  plus_sd_upper_95>0, "TRUE", "FALSE"))

write.csv (selection, snakemake@output[[1]], row.names=F, quote=F)
