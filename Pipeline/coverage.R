library (dplyr)
coverage.df <- read.table(snakemake@input[[1]],stringsAsFactors=F,comment.char = '#') 


coverage.df <- rename (coverage.df, Seg = V1, Pos =V2, Coverage= V3)


filename <-snakemake@input[[1]]
filename2 <- sub(".coverage.csv", "", filename)

filename_vec <- strsplit(filename2, split = "/")[[1]]
coverage.df$sample <- filename_vec[4]


write.table (coverage.df, snakemake@output[[1]], quote=F, row.names=F)

Avg <- mean (coverage.df$Coverage)

summary.df <-  data.frame ("sample"= filename_vec[4], "mean"= Avg)

write.table (summary.df, snakemake@output[[2]], quote=F, row.names=F)
