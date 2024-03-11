library(tidyverse)
library(broom)
library (ggplot2)
library(grid)
library(gridExtra)

snv_corrected <- read.csv ("Results/SNV_with_corrected_frequency_all_mutations.csv", colClasses = c(household = "character", individual= "character", specimen_number= "character"))
meta <- read.csv ("metadata/Vanderbilt_metadata_all_years.csv",  colClasses = c(household = "character", individual= "character", specimen_number= "character"))

# Import information on the number of available sites of each mutation type.
AvailableSites <- read.csv("References/AvailableSitesByType.csv",stringsAsFactors = FALSE)

# remove high snv samples (>50)
high_snv <- read.csv ("Results/high_snv_samples.csv")
meta <- filter (meta, !( CaseID %in% high_snv$Study_ID))

# Import information on coding sequence lengths.
CodingSequenceLengths <- read.csv("References/CodingSequenceLengths.csv", stringsAsFactors = FALSE)


## divergence rate for all samples 
Genome_length <- CodingSequenceLengths %>% 
  filter (! (REGION =="NS") , ! (REGION =="M")) %>%
  group_by (reference) %>% summarise (Length = sum (Length))

DataWithinDiv<- snv_corrected  %>% 
  filter (! (REGION =="NS") , ! (REGION =="M")) %>%
  filter (!( sample %in% high_snv$Study_ID)) %>%
  group_by(reference, sample, mutation_type, PCR_RESULT_1) %>%
  summarize(Div=sum(avg_freq_initial_ref))
# Identify samples that are not represented in this dataframe
# because they have zero variants.
# Add them to the dataframe in the form of a dummy value.

sample_no_ISNV <- meta %>% ungroup () %>% 
  rename (sample = Study_ID) %>%  
  filter(!(sample %in% DataWithinDiv$sample))%>%
  dplyr::select(reference, sample) %>%
  mutate( mutation_type="Non",Div=0)

DataWithinDiv <- rbind (DataWithinDiv, sample_no_ISNV )


# Fill in zero values for samples, genes, and mutation classes
# that have zero reported variants using the complete function.
DataWithinDiv <- DataWithinDiv %>%
  group_by(reference) %>%
  complete(sample, mutation_type, fill=list(Div=0))

# Add metadata about the length of each coding sequence.
DataWithinDiv<- left_join(DataWithinDiv, Genome_length,
                                by=c("reference"))

# Add metadata about the proportion of sites of each mutation type.
DataWithinDiv <- left_join(DataWithinDiv, AvailableSites, 
                                 by=c("mutation_type"))


# Normalize sample divergence based on the number of available sites.
DataWithinDiv <- DataWithinDiv%>%
  mutate(DivPerSite=Div/(Length*PercentSites))

# Normalize sample divergence based on days post infection
meta_sub <- meta %>% rename (sample = Study_ID) %>% select ( sample, Days_post_symp)
DataWithinDiv <- left_join(DataWithinDiv, meta_sub, by = "sample" )  %>%
  mutate (DPI = Days_post_symp +2 )

DataWithinDiv <- DataWithinDiv %>% mutate (DivPerDay =  DivPerSite/DPI)

write.csv(DataWithinDiv, "Results/SampleDivPerSiteGenomeAllSamples.csv",
          quote=FALSE, row.names=FALSE)




## Use lowest ct value if multiple samples per person
meta_ct <- meta %>% group_by (hhsubid) %>% slice_min (subtype_ct)

DataWithinDivCt<- snv_corrected  %>% filter (sample %in% meta_ct$CaseID) %>%
  filter (! (REGION =="NS") , ! (REGION =="M")) %>%
  filter (!( sample %in% high_snv$Study_ID)) %>%
  group_by(reference, sample, mutation_type, PCR_RESULT_1) %>%
  summarize(Div=sum(avg_freq_initial_ref))

sample_no_ISNV_Ct <- meta_ct %>% ungroup () %>% 
  rename (sample = Study_ID) %>%  
  filter(!(sample %in% DataWithinDivCt$sample))%>%
  dplyr::select(reference, sample) %>%
  mutate( mutation_type="Non",Div=0)


DataWithinDivCt <- rbind (DataWithinDivCt, sample_no_ISNV_Ct)

DataWithinDivCt <- DataWithinDivCt %>%
  group_by(reference) %>%
  complete(sample, mutation_type, fill=list(Div=0))

# Add metadata about the length of each coding sequence.
DataWithinDivCt<- left_join(DataWithinDivCt, Genome_length,
                          by=c("reference"))

# Add metadata about the proportion of sites of each mutation type.
DataWithinDivCt <- left_join(DataWithinDivCt, AvailableSites, 
                           by=c("mutation_type"))


# Normalize sample divergence based on the number of available sites.
DataWithinDivCt <- DataWithinDivCt%>%
  mutate(DivPerSite=Div/(Length*PercentSites))

meta_ct <- meta_ct %>% ungroup () %>% rename (sample = CaseID) %>% select (sample, Days_post_symp, cdc_flu_vx, ever_antiviral, age.new, sex, reference, PCR_RESULT_1 )
DataWithinDivCt <- left_join(DataWithinDivCt, meta_ct, by = "sample" )  %>%
  mutate (DPI = Days_post_symp +2 )

DataWithinDivCt<- DataWithinDivCt %>% mutate (DivPerDay =  DivPerSite/DPI)


write.csv(DataWithinDivCt, "Results/SampleDivPerSiteGenomeLowCt.csv",
          quote=FALSE, row.names=FALSE)




#### Divergence by segment###

DataWithinDivSeg <- snv_corrected  %>%
  filter (!( sample %in% high_snv$Study_ID)) %>%
  group_by(reference, sample, GFF_FEATURE, mutation_type) %>%
  summarize(Div=sum(avg_freq_initial_ref))


# Identify samples that are not represented in this dataframe
# because they have zero variants.
# Add them to the dataframe in the form of a dummy value.



sample_no_ISNV <- meta_ct %>% ungroup () %>% 
  filter(!(sample %in% DataWithinDivSeg$sample))%>%
  dplyr::select(reference, sample) %>%
  mutate(GFF_FEATURE="cds-HA", mutation_type="Non",Div=0)

DataWithinDivSeg <- rbind (DataWithinDivSeg, sample_no_ISNV )

# Fill in zero values for samples, genes, and mutation classes
# that have zero reported variants using the complete function.
DataWithinDivSeg <- DataWithinDivSeg %>%
  group_by(reference) %>%
  complete(sample, GFF_FEATURE, mutation_type, fill=list(Div=0))

# Add metadata about the length of each coding sequence.
DataWithinDivSeg<- left_join(DataWithinDivSeg, CodingSequenceLengths,
                          by=c("reference","GFF_FEATURE"))

# Add metadata about the proportion of sites of each mutation type.
DataWithinDivSeg <- left_join(DataWithinDivSeg, AvailableSites, 
                           by=c("mutation_type"))

# Normalize sample divergence based on the number of available sites.
DataWithinDivSeg <- DataWithinDivSeg %>%
  mutate(DivPerSite=Div/(Length*PercentSites))

## add in DPI and remove M/NS
DataWithinDivSeg <- left_join(DataWithinDivSeg, meta_ct, by = "sample" )  %>%
  mutate (DPI = Days_post_symp +2 )

DataWithinDivSeg <- mutate (DataWithinDivSeg, DivPerDay = DivPerSite/DPI) %>%
  filter (! (REGION =="NS") , ! (REGION =="M")) 


write.csv(DataWithinDivSeg, "Results/SampleDivPerSiteSeg.csv",
          quote=FALSE, row.names=FALSE)






## plot by various factors
DataWithinDivCt$cdc_flu_vx <- na_if (DataWithinDivCt$cdc_flu_vx, "Recent vaccination or unknown")

vaccine_sum <- DataWithinDivCt %>% filter (!is.na (DivPerDay)) %>%
  filter (!is.na (cdc_flu_vx)) %>%
  group_by(cdc_flu_vx, mutation_type)%>% 
  summarize (mean = mean (DivPerDay), n= n(), sd =sd (DivPerDay), SE = sd/sqrt(n))


vaccine_mean_plot <- ggplot (vaccine_sum, aes (cdc_flu_vx, mean) )+ 
  geom_pointrange (aes (ymin = mean-SE, ymax = mean+SE, color= mutation_type), 
                   position = position_dodge(width =0.4), show.legend = FALSE) +
  scale_color_manual(values=c("#77427F", "#21A386"))+
  theme_cowplot  (12 )+
  xlab ("Vaccination") + ylab ("")


## anti viral

antiviral_sum <- DataWithinDivCt %>% filter (!is.na (DivPerDay)) %>%
  filter (!is.na (ever_antiviral)) %>%
  group_by(ever_antiviral, mutation_type)%>% 
  summarize (mean = mean (DivPerDay), n= n(), sd =sd (DivPerDay), SE = sd/sqrt(n))

antiviral_sum$ever_antiviral <- as.factor (antiviral_sum$ever_antiviral)
antiviral_mean_plot <- ggplot (antiviral_sum, aes (ever_antiviral, mean) )+ 
  geom_pointrange (aes (ymin = mean-SE, ymax = mean+SE, color= mutation_type), 
                   position = position_dodge(width =0.4), show.legend = FALSE) +
  scale_color_manual(values=c("#77427F", "#21A386"))+
  scale_x_discrete(labels=c("0"="No", "1"= "Yes"))+
  theme_cowplot  (12 )+
  xlab ("Antiviral usage") + ylab ("")

### Age 

## order factors  by age

DataWithinDivCt <- DataWithinDivCt %>% mutate (age_class = ifelse (age.new <= 5, "Child", ifelse 
                                                                   (age.new > 18, "Adult", "Adolescent")))

DataWithinDivCt$age_class <- ordered(DataWithinDivCt$age_class, levels = c("Child", "Adolescent", "Adult"))


### mean +-SE
age_sum <- DataWithinDivCt %>% filter (!is.na (DivPerDay)) %>%
  group_by(age_class, mutation_type)%>% 
  summarize (mean = mean (DivPerDay), n= n(), sd =sd (DivPerDay), SE = sd/sqrt(n))

age_class_mean_plot <- ggplot (age_sum, aes (age_class, mean) )+ 
  geom_pointrange (aes (ymin = mean-SE, ymax = mean+SE, color= mutation_type), 
                   position = position_dodge(width =0.4), show.legend = FALSE) +
  scale_color_manual(values=c("#77427F", "#21A386"))+
  theme_cowplot  (12 )+
  xlab ("Age") + ylab ("")

## subtype 


subtype_sum <- DataWithinDivCt %>% filter (!is.na (DivPerDay)) %>%
  filter (!is.na (PCR_RESULT_1.x)) %>%
  group_by(PCR_RESULT_1.x, mutation_type)%>% 
  summarize (mean = mean (DivPerDay), n= n(), sd =sd (DivPerDay), SE = sd/sqrt(n))

subtype_mean_plot <- ggplot (subtype_sum, aes (PCR_RESULT_1.x, mean) )+ 
  geom_pointrange (aes (ymin = mean-SE, ymax = mean+SE, color= mutation_type), 
                   position = position_dodge(width =0.4), show.legend = FALSE) +
  scale_color_manual(values=c("#77427F", "#21A386"))+
  scale_x_discrete(labels=c("0"="No", "1"= "Yes"))+
  theme_cowplot  (12 )+
  xlab ("Subtype") + ylab ("")


### segment

segment_rate_plot <-  ggplot (  DataWithinDivSeg, aes (REGION, DivPerDay, fill = mutation_type)) + 
  geom_boxplot ( aes (REGION), show.legend = FALSE) +
  scale_fill_manual(values=c("#77427F", "#21A386"))+
  scale_x_discrete(labels=c("NA_"="NA"))+
  theme_cowplot  (12 )+
  xlab ("Age") + ylab ("")



segment_sum <-  DataWithinDivSeg %>% filter (!is.na (DivPerDay)) %>% 
  group_by(REGION, mutation_type)%>% 
  summarize (mean = mean (DivPerDay), n= n(), sd =sd (DivPerDay), SE = sd/sqrt(n))

segment_sum$REGION <- ordered(segment_sum$REGION, levels = c("PB2", "PB1", "PA", "HA", "NP", "NA_"))

segment_rate_mean_plot <-  ggplot (  segment_sum, aes (REGION, mean, fill = mutation_type)) + 
  geom_pointrange (aes (ymin = mean-SE, ymax = mean+SE, color= mutation_type), position = position_dodge(width =0.4), show.legend = FALSE) +
  scale_color_manual(values=c("#77427F", "#21A386"))+
  scale_x_discrete(labels=c("NA_"="NA"))+
  theme_cowplot  (12 )+
  xlab ("Segment") + ylab ("")


### put all figures on one panel 
Upper_plot <- plot_grid (segment_rate_mean_plot, subtype_mean_plot, rel_widths = c(0.75, .5))
lower_plot <-  plot_grid ( antiviral_mean_plot, vaccine_mean_plot, age_class_mean_plot,nrow=1,  rel_widths = c(1, 1, 1.25))
subtype_mean_plot

rate_plot <- plot_grid (Upper_plot, lower_plot, nrow=2 )

rate_plot <- plot_grid (segment_rate_mean_plot, antiviral_mean_plot, age_class_mean_plot, vaccine_mean_plot, nrow=2)
y.grob <- textGrob("Divergence per site per day", 
                   gp=gpar( fontsize=12), rot=90)
grid.arrange(arrangeGrob(rate_plot, left = y.grob))



### stats 

DataWithinDivCtNon <- filter (DataWithinDivCt, mutation_type == "Non")
DataWithinDivCtSyn <- filter (DataWithinDivCt, mutation_type == "Syn")

wilcox.test(DivPerDay ~ cdc_flu_vx, data=DataWithinDivCt ) 
wilcox.test(DivPerDay ~ cdc_flu_vx, data=DataWithinDivCtNon ) 
wilcox.test(DivPerDay ~ cdc_flu_vx, data=DataWithinDivCtSyn ) 


wilcox.test(DivPerDay ~ ever_antiviral, data=DataWithinDivCt ) 
wilcox.test(DivPerDay ~ ever_antiviral, data=DataWithinDivCtNon ) 
wilcox.test(DivPerDay ~ ever_antiviral, data=DataWithinDivCtSyn ) 


wilcox.test(DivPerDay ~ PCR_RESULT_1.x, data=DataWithinDivCt ) 
wilcox.test(DivPerDay ~ PCR_RESULT_1.x, data=DataWithinDivCtNon ) 
wilcox.test(DivPerDay ~ PCR_RESULT_1.x, data=DataWithinDivCtSyn ) 


DataWithinDivCt<- DataWithinDivCt %>% mutate (age_class = ifelse (age.new <=5, "Child", ifelse (age.new >18 , "Adult", "Adolescent")))
DataWithinDivCtNon<- DataWithinDivCtNon %>% mutate (age_class = ifelse (age.new <=5, "Child", ifelse (age.new >18 , "Adult", "Adolescent")))
DataWithinDivCtSyn<- DataWithinDivCtSyn %>% mutate (age_class = ifelse (age.new <=5, "Child", ifelse (age.new >18 , "Adult", "Adolescent")))


kruskal.test(DivPerDay ~ age_class, data = DataWithinDivCt)
kruskal.test(DivPerDay ~ age_class, data = DataWithinDivCtNon)
kruskal.test(DivPerDay ~ age_class, data = DataWithinDivCtSyn)

DataWithinDivSegNon <- filter (DataWithinDivSeg, mutation_type == "Non")
DataWithinDivSegSyn<- filter (DataWithinDivSeg, mutation_type == "Syn")


kruskal.test(DivPerDay ~ REGION, data = DataWithinDivSeg)
kruskal.test(DivPerDay ~ REGION, data = DataWithinDivSegNon)
