library (plyr)
library (dplyr)
library (stringr)
library (Biostrings)

variant.df <- read.table(snakemake@input[[1]],stringsAsFactors=F,comment.char = '#', header=T)
names(variant.df)[10:29] <- c('REF_DP_1', 'REF_RV_1', 'REF_QUAL_1', 'ALT_DP_1', 'ALT_RV_1', 'ALT_QUAL_1', 'ALT_FREQ_1', 'TOTAL_DP_1', 'PVAL_1', 'PASS_1','REF_DP_2', 'REF_RV_2', 'REF_QUAL_2', 'ALT_DP_2', 'ALT_RV_2', 'ALT_QUAL_2', 'ALT_FREQ_2', 'TOTAL_DP_2', 'PVAL_2', 'PASS_2')

##basic filtering
variant.df<-filter (variant.df, !str_detect (ALT, "[+-]")) # remove indels 
variant.df<-mutate(variant.df,mutation=paste0(REGION,"_",REF,POS,ALT))
#variant.df<-subset(variant.df, PASS_1=="TRUE" & PASS_2=="TRUE") # only keep variants that pass p value (0.05)
variant.df <- filter (variant.df, PVAL_1 <=1e-5 &  PVAL_2  <=1e-5) # only keep variants that have a p value  less than 10^-5 
variant.df <- filter (variant.df, TOTAL_DP_1 >=400, TOTAL_DP_2 >=400) # total depth is  equal to or greater than 400. 

reference.fasta<-snakemake@input[[2]] 
segments <- fasta.seqlengths(reference.fasta)
regions.bed <- data.frame(chr = gsub("[ ].*","", names(segments)), start=12, stop=segments-13, row.names=NULL) # the univeral primers are 12 and 13 bp long
regions.bed<-mutate(regions.bed,chr=as.character(chr))

primer.cut<-function(x){ # a helper function to do this
  chr<-unique(x$REGION)
  start<-regions.bed$start[match(x$REGION,regions.bed$chr)]
  stop<-regions.bed$stop[match(x$REGION,regions.bed$chr)]
  subset(x,POS>start & POS<stop)
}

variant.df<-ddply(variant.df,~REGION,primer.cut) # removes any variants within  priming sites 


## remove variants where there isn't 0.25% frequency in both and 0.5% frequency in at least one sample
 
#variant.df <- filter (variant.df,ALT_FREQ_1 >=0.0025 & ALT_FREQ_2 >=0.0025, ALT_FREQ_1 <=0.9975 & ALT_FREQ_2 <=0.9975)
variant.df <- filter (variant.df, ALT_FREQ_1 >=0.005 & ALT_FREQ_2 >=0.005)
variant.df <- filter (variant.df, ALT_FREQ_1 <=0.995 & ALT_FREQ_2 <=0.995)

filename <-snakemake@input[[1]]
filename2 <- sub(".merged.tsv", "", filename)

filename_vec <- strsplit(filename2, split = "/")[[1]]
if (dim(variant.df)[1] != 0) {
variant.df$sample <- filename_vec[3]
}

write.table (variant.df, snakemake@output[[1]], quote=F, row.names=F)
