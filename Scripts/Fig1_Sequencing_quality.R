library (tidyr)
library (ggplot2)

#read in per pos coverage
cov_HongKong <- read.table ("Processed_data/HongKong14_A/Coverage.all", header=T)
cov_Michigan <- read.table ("Processed_data/Michigan15_A/Coverage.all", header =T)
cov_Brisbane <- read.table ("Processed_data/Brisbane18_A/Coverage.all", header=T)
cov_Singapore <- read.table ("Processed_data/Singapore16_A/Coverage.all", header=T)

#read in length to add  and join with coverage
length_HongKong <- read.csv ("metadata/HongKong_length_to_add.csv")
length_Michigan <- read.csv ("metadata/Michigan_length_to_add.csv")
length_Brisbane <- read.csv ("metadata/Brisbane_length_to_add.csv")
length_Singapore <- read.csv ("metadata/Singapore_length_to_add.csv")


cov_HongKong<- full_join (cov_HongKong, length_HongKong, by = "Seg")
cov_Michigan <- full_join (cov_Michigan, length_Michigan, by = "Seg")
cov_Brisbane <- full_join (cov_Brisbane, length_Brisbane, by = "Seg")
cov_Singapore <- full_join (cov_Singapore, length_Singapore, by = "Seg")

cov <- bind_rows ( cov_HongKong, cov_Michigan, cov_Brisbane,  cov_Singapore)
cov <- mutate (cov, concat.pos = Pos+ Length_to_add)
# Create concate pos

###### Functions #####
#coverage



slide <- function(cov.df, setup.df)
{
  Coverage = rep(NA, nrow(setup.df))
  for(i in 1:nrow(setup.df))
  {
    s = setup.df$starts[i]
    e = setup.df$ends[i]
    subset(cov.df, concat.pos >= s & concat.pos < e, select = c(Coverage)) -> position
    mean(position$Coverage) -> Coverage[i]
  }
  out <- data.frame(mean = Coverage, concat.pos = setup.df$concat.pos, Seg = setup.df$Seg)
  out$sample = unique(cov.df$sample)
  return(out)
}

cov_plot <- function(cov.df, title)
{
  cov.df %>% group_by(Seg) %>% summarize(first = min(concat.pos), last = max(concat.pos)) %>% plyr::adply(1,function(x) data.frame(starts = seq(x$first,x$last,by=300))) %>% mutate(ends = ifelse(starts + 300 < last, starts + 300, last)) -> setup
  setup %>% select(starts, ends) -> setup_means
  setup$concat.pos <- apply(setup_means, 1, function(x) mean(x))
  plyr::ddply(cov.df, ~sample, slide, setup) -> cov.slid.df
  
  x.labels <- plyr::ddply(cov.slid.df,~Seg,plyr::summarize,concat.pos=concat.pos[which(abs(concat.pos-mean(concat.pos))==(min(abs(concat.pos-mean(concat.pos)))))])
  x.labels <- plyr::ddply(x.labels, ~Seg, function(x) return(x[1,]))
  
  cov.plot <- ggplot(cov.slid.df, mapping = aes(x = as.factor(concat.pos), y = mean)) + geom_boxplot(fill="white")
  cov.plot <- cov.plot + ggtitle(title) + ylab("Read depth") + scale_x_discrete(labels = x.labels$Seg, breaks = x.labels$concat.pos) + xlab("Genome Position")
  cov.plot <- cov.plot + theme(axis.title.y = element_text(vjust=1.2))
  cov.plot <- cov.plot + theme(legend.position = "none") + theme_classic()
  cov.plot <- cov.plot + theme(text = element_text(size = 12), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
  return(cov.plot)
}
# ======================== Plot ===================
cov_plot(cov, title = "") -> coverage.plot.all



## Allele frequency between replicates including an insert of 0 to 0.1 frequency
snv<- read.csv ("Processed_data/Secondary_processing/SNV_with_meta_data.csv",  colClasses = c(household = "character", individual= "character", specimen_number= "character"))


rep_allele_freq_plot <- ggplot (snv, aes (ALT_FREQ_1 , ALT_FREQ_2 )) + geom_point (alpha =.5) + 
  xlab (" Allele Frequency in Replicate 1")+ ylab ("Allele Frequency in Replicate 2") +
  theme_cowplot (12)

snv_rep_plot <- ggplot (snv, aes (ALT_FREQ_1, ALT_FREQ_2))+ geom_point(alpha =0.5, size =0.8, shape=16)+ 
  xlab ("iSNV Frequency in Replicate 1")+ ylab ("iSNV Frequency in Replicate 2")+
  theme_bw()+ xlim (0, 1)+ ylim (0,1) 

plot_insert <- ggplot (snv, aes (ALT_FREQ_1, ALT_FREQ_2))+ geom_point(size=0.3, alpha =0.5, shape=16)+ 
  theme_cowplot(12) + theme (panel.background = element_rect(fill = "white"))+
  scale_x_continuous(breaks = c (0, 0.05, 0.1), limits = c(0, 0.1)) + 
  scale_y_continuous(breaks= c (0, 0.05, 0.1),limits = c(0, 0.1)) + 
  xlab ("")+ ylab ("")

snv_rep_plot_inset <- snv_rep_plot + inset_element(plot_insert, 0.36, .36, .97, .99)



#### Codon Position

snv_meta_corrected <- read.csv ("Results/SNV_with_corrected_frequency_all_mutations.csv",  colClasses = c(household = "character", individual= "character", specimen_number= "character"))

## brisbane - get codon position for each iSNV based on gff file
brisbane_Snv <- filter (snv_meta_corrected, reference== "Brisbane")
brisbane_Snv <- mutate (brisbane_Snv, peptide_position = ifelse (REGION== "PB2", POS- 27, ifelse(
  REGION== "PB1", POS-24, ifelse (
    REGION == "PA", POS-24, ifelse (
      REGION == "HA", POS-32, ifelse (
        REGION == "NP", POS- 45, ifelse (
          REGION == "NA_", POS-20, "NA"
        )))))))
brisbane_Snv$peptide_position <- as.numeric (brisbane_Snv$peptide_position)
brisbane_Snv <- mutate (brisbane_Snv,codon_position = ifelse (peptide_position%%3==0 ,3, ifelse (
  (peptide_position+1)%%3 == 0, 2, ifelse(
    (peptide_position+2)%%3 == 0, 1, NA)
)))

## singapore 

Singapore_Snv <- filter (snv_meta_corrected, reference== "Singapore")
Singapore_Snv <- mutate (Singapore_Snv, peptide_position = ifelse (REGION== "PB2", POS- 27, ifelse(
  REGION== "PB1", POS-24, ifelse (
    REGION == "PA", POS-24, ifelse (
      REGION == "HA", POS-29, ifelse (
        REGION == "NP", POS- 46, ifelse (
          REGION == "NA_", POS-19, "NA"
        )))))))
Singapore_Snv$peptide_position <- as.numeric (Singapore_Snv$peptide_position)
Singapore_Snv <- mutate (Singapore_Snv,codon_position = ifelse (peptide_position%%3==0 ,3, ifelse (
  (peptide_position+1)%%3 == 0, 2, ifelse(
    (peptide_position+2)%%3 == 0, 1, NA)
)))

### HongKong

HongKong_Snv <- filter (snv_meta_corrected, reference== "HongKong")
HongKong_Snv <- mutate (HongKong_Snv, peptide_position = ifelse (REGION== "PB2", POS- 27, ifelse(
  REGION== "PB1", POS-24, ifelse (
    REGION == "PA", POS-24, ifelse (
      REGION == "HA", POS-29, ifelse (
        REGION == "NP", POS- 45, ifelse (
          REGION == "NA_", POS-19, "NA"
        )))))))
HongKong_Snv$peptide_position <- as.numeric (HongKong_Snv$peptide_position)
HongKong_Snv <- mutate (HongKong_Snv,codon_position = ifelse (peptide_position%%3==0 ,3, ifelse (
  (peptide_position+1)%%3 == 0, 2, ifelse(
    (peptide_position+2)%%3 == 0, 1, NA)
)))


## Michigan
Michigan_Snv <- filter (snv_meta_corrected, reference== "Michigan")
Michigan_Snv <- mutate (Michigan_Snv, peptide_position = ifelse (REGION== "PB2", POS- 27, ifelse(
  REGION== "PB1", POS-24, ifelse (
    REGION == "PA", POS-24, ifelse (
      REGION == "HA", POS-32, ifelse (
        REGION == "NP", POS- 45, ifelse (
          REGION == "NA_", POS-20, "NA"
        )))))))
Michigan_Snv$peptide_position <- as.numeric (Michigan_Snv$peptide_position)
Michigan_Snv <- mutate (Michigan_Snv,codon_position = ifelse (peptide_position%%3==0 ,3, ifelse (
  (peptide_position+1)%%3 == 0, 2, ifelse(
    (peptide_position+2)%%3 == 0, 1, NA)
)))

## combine reference specific 
SNV_binned <- bind_rows (brisbane_Snv, HongKong_Snv)
SNV_binned <- bind_rows (SNV_binned, Singapore_Snv)
SNV_binned <- bind_rows (SNV_binned, Michigan_Snv)

## round frequency to nearest 0.05
library (plyr)
SNV_binned  <- SNV_binned  %>% mutate (avg_freq_bin = round_any (avg_freq_initial_ref, 0.05))
detach("package:plyr", unload = TRUE)

## get proportion of codon positions for each frequency bin
SNV_binned_sum <- SNV_binned  %>% count  (avg_freq_bin) %>%  rename("total_mutation" = "n")
SNV_binned_sum_1 <- SNV_binned %>% count (codon_position, avg_freq_bin) %>% filter (codon_position == 1) %>% rename ("Codon_count" = "n") 
SNV_binned_sum_2 <- SNV_binned %>% count (codon_position, avg_freq_bin) %>% filter (codon_position == 2) %>% rename ("Codon_count" = "n") 
SNV_binned_sum_3 <- SNV_binned %>% count (codon_position, avg_freq_bin) %>% filter (codon_position == 3) %>% rename ("Codon_count" = "n") 


SNV_binned_sum_1 <- left_join (SNV_binned_sum  , SNV_binned_sum_1 , by = c( "avg_freq_bin"))

SNV_binned_sum_2  <- left_join (SNV_binned_sum  , SNV_binned_sum_2 , by = c( "avg_freq_bin"))
SNV_binned_sum_3 <- left_join (SNV_binned_sum  , SNV_binned_sum_3 , by = c("avg_freq_bin"))

SNV_binned_sum <- bind_rows(SNV_binned_sum_1, SNV_binned_sum_2, SNV_binned_sum_3)

SNV_binned_sum  <- filter (SNV_binned_sum , !is.na(codon_position)) %>% replace (is.na (.),0) %>% mutate (prop_codon = Codon_count/total_mutation)


### plot codon position
codon_plot <- ggplot (SNV_binned_sum , aes (avg_freq_bin, prop_codon )) + 
  geom_line ( size = 0.75,  aes (color = as.factor (codon_position))) + 
  xlim (0, 0.5) + theme_cowplot (12) + theme (legend.position="none")+ scale_color_manual(values =c ("#F9766D", "#F2AD00", "#5BBCD6"))+
  xlab ("iSNV Frequency") + ylab ("Proportion")



## combine plots
bottom_plot <- plot_grid(
  snv_rep_plot_inset,  codon_plot, align = "h", axis = "bt", 
  ncol = 2, labels = c('B', 'C'), vjust = 0.5)

plot_grid (coverage.plot.all, NULL, bottom_plot, nrow =3, labels = c('A', '', ''), rel_heights = c(1,.1, 1.2))


