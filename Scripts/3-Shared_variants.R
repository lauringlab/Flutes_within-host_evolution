library (tidyr)
library (ggplot2)
library (dplyr)
library(grid)
library(gridExtra)


## read in data
snv_meta <- read.csv ("Processed_data/Secondary_processing/SNV_with_meta_data.csv", colClasses = c(household = "character", individual= "character", specimen_number= "character"))


## Keep only a single isntance of an iSNV per household, and filter based on number of individuals an iSNV is found in
variants_count <-  snv_meta %>% distinct (hhid, mutation,.keep_all=T) %>%
  count (mutation, reference, concat.pos, mutation_type, REGION, PCR_RESULT_1, POS) %>% 
  filter(case_when(reference=='Brisbane' ~ n >= 4,
                   reference=='Michigan' ~ n >= 2,
                   reference == "HongKong" ~ n >=2,
                   reference == "Singapore" ~ n >=3))
               

variants_count <- unite (variants_count, reference, mutation,col="reference_mutation", sep = "-", remove=F)
shared_variants <- snv_meta %>% unite (reference, mutation,col="reference_mutation", sep = "-", remove=F) %>% filter (reference_mutation %in% variants_count$reference_mutation)
write.csv (shared_variants, "Results/shared_variants.csv", row.names=F, quote=F)


## repeat after high iSNV samples are removed
variants_count_high <-  snv_meta %>% filter (!(sample %in% high_snv_sample$sample)) %>%
  distinct (hhid, mutation,.keep_all=T) %>%
  count (mutation, reference, concat.pos, mutation_type, REGION, PCR_RESULT_1, POS) %>% 
  filter(case_when(reference=='Brisbane' ~ n >= 4,
                   reference=='Michigan' ~ n >= 2,
                   reference == "HongKong" ~ n >=2,
                   reference == "Singapore" ~ n >=3))

variants_count_high <- unite (variants_count_high, reference, mutation,col="reference_mutation", sep = "-", remove=F)
shared_variants_no_high <- snv_meta %>% unite (reference, mutation,col="reference_mutation", sep = "-", remove=F) %>% filter (reference_mutation %in% variants_count_high$reference_mutation)
write.csv (shared_variants_no_high, "Results/shared_variants_wo_high_snv_samples.csv", row.names=F, quote=F)


## antigenic sties
antigenic_sites <- read.csv ("References/HA_antigenic_sites_bp.csv")

shared_H3N2_antigenic <- variants_count %>% filter (PCR_RESULT_1 == "FluAH3", REGION == "HA", POS %in% antigenic_sites$H3N2) 
write.csv (shared_H3N2_antigenic, "Results/H3N2_antigenic_sites_in_3_or_more_individuals.csv", row.names=F, quote=F)

shared_H1N1_antigenic <- variants_count %>% filter (PCR_RESULT_1 == "FluAH1N1", REGION == "HA", POS %in% antigenic_sites$H1N1) # no antigenic sites present

## summarize and plot 
shared_variants_sum <- shared_variants %>% ungroup () %>% group_by (reference, mutation) %>% summarize (mean_freq= mean(avg_freq), 
                        median = median (avg_freq), min = min (avg_freq), max = max (avg_freq), 
                        sd_neg = mean (avg_freq)- sd (avg_freq), sd_pos = mean (avg_freq)+sd (avg_freq))
                        
shared_variants_sum <- inner_join (shared_variants_sum,variants_count , by = c("mutation", "reference"))


shared_variants_sum_Michigan <- filter (shared_variants_sum , reference == "Michigan")
Michigan_shared_variant_plot <- ggplot (shared_variants_sum_Michigan, aes (n, mean_freq, color= mutation_type)) + 
  geom_jitter (alpha=0.65) + 
  scale_y_continuous(trans='log10', limits = c (0.003, 0.2)) + xlim (0,14) +
  xlab ("Number of individuals sharing an iSNV") + ylab ("Mean Frequency")+
  facet_wrap (~factor (REGION, levels= c ("PB2", "PB1", "PA", "NP", "HA", "NA_", "M", "NS")), drop=F, ncol =4) +
  theme_bw () + theme (legend.position = "none" , plot.title = element_text(  face = "bold", size = 10)) +
  scale_color_manual(values=c("#481468", "#21A386")) +
  ggtitle ("H1N1 2017 & 2018")
 

shared_variants_sum_Brisbane <- filter (shared_variants_sum , reference == "Brisbane")
Brisbane_shared_variant_plot <- ggplot (shared_variants_sum_Brisbane, aes (n, mean_freq, color= mutation_type)) + 
   geom_jitter (alpha=0.65) + 
   scale_y_continuous(trans='log10', limits = c (0.003, 0.2)) +
   xlab ("Number of individuals sharing an iSNV") + ylab ("")+
   facet_wrap (~factor (REGION, levels= c ("PB2", "PB1", "PA", "NP", "HA", "NA_", "M", "NS")), drop=F, ncol =4) +
   theme_bw () + theme (legend.position = "none", plot.title = element_text( face = "bold", size =10))+
   scale_color_manual(values=c("#481468", "#21A386"))+
  ggtitle ("H1N1 2019")


shared_variants_sum_HongKong <- filter (shared_variants_sum , reference == "HongKong")
HongKong_shared_variant_plot <- ggplot (shared_variants_sum_HongKong, aes (n, mean_freq, color= mutation_type)) + 
  geom_jitter (alpha=0.65) + 
  scale_y_continuous(trans='log10', limits = c (0.003, 0.2)) + xlim (0,11) +
  xlab ("") + ylab ("Mean Frequency")+
  facet_wrap (~factor (REGION, levels= c ("PB2", "PB1", "PA", "NP", "HA", "NA_", "M", "NS")), drop=F, ncol =4) +
  theme_bw () + theme (legend.position = "none", plot.title = element_text(  face = "bold", size =10)) +
  scale_color_manual(values=c("#481468", "#21A386")) +
  ggtitle ("H3N2 2017")

shared_variants_sum_Singapore <- filter (shared_variants_sum , reference == "Singapore")
Singapore_shared_variant_plot <- ggplot (shared_variants_sum_Singapore, aes (n, mean_freq, color= mutation_type)) + 
  geom_jitter (alpha=0.65) + 
  scale_y_continuous(trans='log10', limits = c (0.003, 0.2)) + xlim (2,15) +
  xlab ("") + ylab ("")+
  facet_wrap(~factor (REGION, levels= c ("PB2", "PB1", "PA", "NP", "HA", "NA_", "M", "NS")), drop=F, ncol =4) +
  theme_bw () +  theme (legend.position = "none", plot.title = element_text( face = "bold", size =10)) +
  scale_color_manual(values=c("#481468", "#21A386")) +
  ggtitle ("H3N2 2018")

 
shared_variants_plot_wrap <- plot_grid (HongKong_shared_variant_plot, Singapore_shared_variant_plot , Michigan_shared_variant_plot, Brisbane_shared_variant_plot, ncol=2 )

ggsave ("Results/Plots/shared_variants_plot_wrap.pdf", plot =shared_variants_plot_wrap, height = 5, width = 8)



## summarize and plot shared variants without high SNV samples

shared_variants_sum_no_high <- shared_variants_no_high %>% group_by (reference, mutation) %>% summarize (mean_freq= mean(avg_freq), 
                                                                                         median = median (avg_freq), min = min (avg_freq), max = max (avg_freq), 
                                                                                         sd_neg = mean (avg_freq)- sd (avg_freq), sd_pos = mean (avg_freq)+sd (avg_freq))

shared_variants_sum_no_high <- inner_join (shared_variants_sum_no_high,variants_count_high , by = c("mutation", "reference"))



shared_variants_sum_Michigan_no_hign <- filter (shared_variants_sum_no_high , reference == "Michigan")
Michigan_shared_variant_plot_no_high <- ggplot (shared_variants_sum_Michigan_no_hign, aes (n, mean_freq, color= mutation_type)) + 
  geom_jitter (alpha=0.65) + 
  scale_y_continuous(trans='log10', limits = c (0.003, 0.2)) + xlim (0,13) +
  xlab ("Number of individuals sharing an iSNV") + ylab ("Mean Frequency")+
  facet_wrap (~factor (REGION, levels= c ("PB2", "PB1", "PA", "NP", "HA", "NA_", "M", "NS")), drop=F, ncol =4) +
  theme_bw () + theme (legend.position = "none" , plot.title = element_text(  face = "bold", size =10)) +
  scale_color_manual(values=c("#481468", "#21A386")) +
  ggtitle ("H1N1 2017 & 2018")


shared_variants_sum_Brisbane_no_high <- filter (shared_variants_sum_no_high , reference == "Brisbane")
Brisbane_shared_variant_plot_no_high <- ggplot (shared_variants_sum_Brisbane_no_high, aes (n, mean_freq, color= mutation_type)) + 
  geom_jitter (alpha=0.65) + 
  scale_y_continuous(trans='log10', limits = c (0.003, 0.2)) +
  xlab ("Number of individuals sharing an iSNV") + ylab ("")+
  facet_wrap (~factor (REGION, levels= c ("PB2", "PB1", "PA", "NP", "HA", "NA_", "M", "NS")), drop=F, ncol =4) +
  theme_bw () + theme (legend.position = "none", plot.title = element_text( face = "bold", size =10))+
  scale_color_manual(values=c("#481468", "#21A386"))+
  ggtitle ("H1N1 2019")


shared_variants_sum_HongKong_no_high <- filter (shared_variants_sum_no_high, reference == "HongKong")
HongKong_shared_variant_plot_no_high <- ggplot (shared_variants_sum_HongKong_no_high, aes (n, mean_freq, color= mutation_type)) + 
  geom_jitter (alpha=0.65) + 
  scale_y_continuous(trans='log10', limits = c (0.003, 0.2)) +
  xlab ("") + ylab ("Mean Frequency")+
  facet_wrap (~factor (REGION, levels= c ("PB2", "PB1", "PA", "NP", "HA", "NA_", "M", "NS")), drop=F, ncol =4) +
  theme_bw () + theme (legend.position = "none", plot.title = element_text(  face = "bold", size =10)) +
  scale_color_manual(values=c("#481468", "#21A386")) +
  ggtitle ("H3N2 2017")

shared_variants_sum_Singapore_no_high <- filter (shared_variants_sum_no_high , reference == "Singapore")
Singapore_shared_variant_plot_no_high <- ggplot (shared_variants_sum_Singapore_no_high, aes (n, mean_freq, color= mutation_type)) + 
  geom_jitter (alpha=0.65) + 
  scale_y_continuous(trans='log10', limits = c (0.003, 0.2)) +  xlim (2,15) +
  xlab ("") + ylab ("")+
  facet_wrap(~factor (REGION, levels= c ("PB2", "PB1", "PA", "NP", "HA", "NA_", "M", "NS")), drop=F, ncol =4) +
  theme_bw () +  theme (legend.position = "none", plot.title = element_text( face = "bold", size =10)) +
  scale_color_manual(values=c("#481468", "#21A386")) +
  ggtitle ("H3N2 2018")


shared_variants_plot_wrap_no_high <- plot_grid (HongKong_shared_variant_plot_no_high, Singapore_shared_variant_plot_no_high , Michigan_shared_variant_plot_no_high, Brisbane_shared_variant_plot_no_high, ncol=2 )

ggsave ("Results/Plots/shared_variants_plot_wrap_no_high.pdf", plot =shared_variants_plot_wrap_no_high, height = 5, width = 8)


##plot NS 2017 H3N2 - look at haplotype

HongKong_NS <- filter (shared_variants_no_high, reference == "HongKong", REGION == "NS")
ggplot (HongKong_NS, aes (POS, as.factor (sample)))+ geom_point () + 
  theme_bw () + xlab ("Position in NS") + ylab ("Sample")



ggplot (HongKong_NS, aes (POS, as.factor (sample)))+ geom_point (size=1)+
  +   theme_bw () + xlab ("Position in NS") + ylab ("Sample") + xlim (0,890)
