library (dplyr)
library (ggplot2)
library (cowplot)
library (tidyr)

# read in selection coefficients
H1N1_selection_Coef <- read.csv ("Results/WFABC/FluAH1N1/H1N1_selection_coefficients.csv")
H3N2_selection_Coef <- read.csv ("Results/WFABC/FluAH3/H3N2_selection_coefficients.csv")

# read in other mutational data 
snv <- read.csv ("Results/SNV_with_corrected_frequency.csv", colClasses = c(household = "character", individual= "character", specimen_number= "character"))
snv_ps <- select (snv, hhsubid, hhsubid_pos, REGION, POS, mutation_type, initial_mutation, concat.pos) %>% distinct (hhsubid_pos, .keep_all=T)

meta <- read.csv ("metadata/Vanderbilt_metadata_all_years.csv", colClasses = c(household = "character", individual= "character", specimen_number= "character"))
meta_sub <- select (meta, hhsubid, cdc_flu_vx, ever_antiviral, age.new, sex, age_class)
meta_sub <- distinct (meta_sub, hhsubid, .keep_all=T)



#filter for positive selection coefficients- 95 posterior density cannot cross 0 for any of the 3 sample sizes
H1N1_pos <-filter (H1N1_selection_Coef, mean_pos_selection == "TRUE",
                   minus_sd_pos_selection == "TRUE", plus_sd_pos_selection == "TRUE")
#H1N1_mean_ne_pos <- filter (H1N1_selection_Coef, mean_pos_selection == "TRUE")
H1N1_pos <- mutate (H1N1_pos, subtype= "FluAH1N1")
H1N1_pos <-left_join (H1N1_pos, snv_ps, by = "hhsubid_pos") 
write.csv (H1N1_pos, "Results/H1N1_pos_selection.csv", row.names=F, quote=F)

H3N2_pos <-filter (H3N2_selection_Coef, mean_pos_selection == "TRUE",
                   minus_sd_pos_selection == "TRUE", plus_sd_pos_selection == "TRUE")
#H3N2_mean_ne_pos <- filter (H3N2_selection_Coef, mean_pos_selection == "TRUE") 2 positions where 95% posterior crossed zero for the 2 other Nes but not for the mean Ne
H3N2_pos <- mutate (H3N2_pos, subtype= "FluAH3")
H3N2_pos <-left_join (H3N2_pos, snv_ps, by = "hhsubid_pos")
write.csv (H3N2_pos, "Results/H3N2_pos_selection.csv", row.names=F, quote=F)

#combine subtypes and plot selection coefficients
pos_selection_coef <- bind_rows (H3N2_pos, H1N1_pos)

## check for antigenic sites
antigenic_sites <- read.csv ("References/HA_antigenic_sites_bp.csv")

H3N2_antigenic <- filter (H3N2_pos, REGION == "HA", POS %in% antigenic_sites$H3N2)
write.csv (H3N2_antigenic, "Results/H3N2_Positive_selection_antigenic.csv", row.names=F, quote=F)

H1N1_antigenic <- filter (H1N1_pos, REGION == "HA", POS %in% antigenic_sites$H1N1)
write.csv (H1N1_antigenic, "Results/H1N1_Positive_selection_antigenic.csv", row.names=F, quote=F)

## see if positively selected alleles became major allele or were found in mutliple individuals
pos_selection_coef <- bind_rows (H3N2_pos, H1N1_pos)
major_minor_allele <- read.csv ("Results/major_minor_allele_trajectory.csv", 
                                colClasses = c(household = "character", individual= "character", specimen_number= "character"))
major_minor_allele_pos_selection <- filter  (pos_selection_coef, hhsubid_pos %in%  major_minor_allele$hhsubid_pos )

write.csv (major_minor_allele_pos_selection, "Results/major_minor_allele_pos_selection.csv", row.names=F, quote =F)


shared_variants <- read.csv ("Results/shared_variants.csv")
shared_variants_H1N1 <- filter (shared_variants, PCR_RESULT_1 == "FluAH1N1")
H1N1_pos_shared <- filter (H1N1_pos, initial_mutation %in% shared_variants_H1N1$mutation) # none shared

shared_variants_H3N2 <- filter (shared_variants, PCR_RESULT_1 == "FluAH3")
H3N2_pos_shared <- filter (H3N2_pos, initial_mutation %in% shared_variants_H3N2$mutation) # 1 shared
write.csv (H3N2_pos_shared, "Results/shared_variants_w_pos_S.csv", row.names=F, quote=F)


shared_major_minor_allele_h1n1 <- filter (shared_variants_H1N1, mutation %in% major_minor_allele$mutation) # found 11
write.csv (shared_major_minor_allele_h1n1, "Results/H1N1_shared_major_minor_allele.csv", row.names=F, quote=F)

shared_major_minor_allele_h3n2 <- filter (shared_variants_H3N2, mutation %in% major_minor_allele$mutation) # found 7
write.csv (shared_major_minor_allele_h3n2, "Results/H3N2_shared_major_minor_allele.csv", row.names=F, quote=F)


###Plot

H1N1_plot <- ggplot (H1N1_pos, aes (concat.pos, mean_s))+ geom_point(aes (color=mutation_type))+ 
  geom_segment( aes(x=concat.pos, xend=concat.pos, y=0, yend=mean_s, color = mutation_type)) +
  scale_color_manual(values=c("#481468", "#21A386")) + 
  scale_x_continuous (breaks = c (1170, 3511, 5798, 7803, 9474, 10986, 12229, 13181), labels = c ("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"), lim = c (0, 13631)) +
  ggtitle ("H1N1") + ylab ("Selection Coefficient") + xlab ("") + 
  theme_cowplot (12) + theme (plot.title = element_text(size=10), axis.text=element_text(size=10), axis.title=element_text(size=10) ) 


H3N2_plot <- ggplot (H3N2_pos, aes (concat.pos, mean_s))+ geom_point(aes (color=mutation_type))+ 
  geom_segment( aes(x=concat.pos, xend=concat.pos, y=0, yend=mean_s, color = mutation_type)) +
  scale_color_manual(values=c("#481468", "#21A386")) + 
  scale_x_continuous (breaks = c (1170, 3511, 5798, 7796, 9461, 10978, 12224, 13182), labels = c ("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"), lim = c (0, 13631)) +
   ggtitle ("H3N2") + ylab ("Selection Coefficient") + xlab ("Genome Position") + 
  theme_cowplot (12) + theme (plot.title = element_text(size=10), axis.text=element_text(size=10), axis.title=element_text(size=10) ) 

## allele trajectory plot 


Major_minor <- read.csv ("Results/Major_minor_allele_trajectory.csv", colClasses = c(household = "character", individual= "character", specimen_number= "character")) 

## add in dummy color variable
color_var <-read.csv ("metadata/color_var_allele_traj.csv")
Major_minor <- left_join (Major_minor,color_var, by = c("hhsubid", "PCR_RESULT_1"))

# split by subtype
Major_minor$PCR_RESULT_1 <- factor(Major_minor$PCR_RESULT_1, levels = c("FluAH1N1", "FluAH3"), 
                                   labels = c("H1N1", "H3N2"))

Major_minor_plot <- ggplot (Major_minor, aes (Days_post_symp_hh, avg_freq_initial_ref, group =hhsubid_pos )) + 
  geom_line (aes ( color = as.factor (color_var),linetype = mutation_type)) + 
  xlab ("Days post symptom onset") + ylab ("iSNV frequency")+ 
  theme_bw () +  theme (legend.position="none")+ scale_x_continuous(breaks=seq(0,9,1)) +
  facet_grid (rows=vars (PCR_RESULT_1) )


## arrange plots

selection_coefficient_plot <- plot_grid(
  H1N1_plot  + theme(legend.position="none"),
  H3N2_plot + theme(legend.position="none"),
  hjust = -1,
  nrow = 2
)

allele_traj_plot <- plot_grid (Major_minor_plot, selection_coefficient_plot, nrow=2, labels = c('A', 'B'),label_size = 12)








