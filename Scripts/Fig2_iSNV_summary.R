## overview analyses
library (dplyr)
library (ggplot2)
library (cowplot)


## read in data 
snv_meta <- read.csv ("Processed_data/Secondary_processing/SNV_with_meta_data.csv", colClasses = c(household = "character", individual= "character", specimen_number= "character"))


# Number of samples per individual
samples_per_individual <- snv_meta %>% distinct (Study_ID, .keep_all=T) %>% count (hhsubid, reference, PCR_RESULT_1)
ggplot (samples_per_individual, aes (as.factor (n)))+ geom_bar () + xlab ("Number of samples per individual")+
  ylab ("Number of individuals") +theme_cowplot (12) +
  theme(axis.text=element_text(size=10))
ggsave ("Results/Plots/Num_samples_per_individual.pdf")



### Figure 2 here
## count the number of iSNV per sample
snv_meta_count <- count (snv_meta, Study_ID,PCR_RESULT_1, cdc_flu_vx, ever_antiviral, age.new, sex, hhsubid, hhid, Days_post_symp, age_class)
sample_no_isnv <- filter (meta, ! (Study_ID %in% snv_meta$sample))

snv_meta_count [nrow(snv_meta_count) + 1,] = c( "1908901505", "FluAH1N1", "Vaccinated","1","42", "1", "1908901", "19089", "9", "Adult", "0" )
snv_meta_count[nrow(snv_meta_count) + 1,] = c( "1909901502", "FluAH1N1", "Vaccinated","1", "4", "1", "1909901", "19099","7","Young_child", "0" )
snv_meta_count [nrow(snv_meta_count) + 1,] = c( "1807502502", "FluAH3",  "Vaccinated","1", "4", "0", "1807502", "18075",NA, "Young_child", "0" )
snv_meta_count [nrow(snv_meta_count) + 1,] = c( "1911502501", "FluAH1N1", "Vaccinated","1", "11", "1", "1911502","19115","1", "Adolescent", "0" )


snv_meta_count$n <-as.numeric (snv_meta_count$n)
snv_per_sample_plot <-ggplot (snv_meta_count, aes (n)) + geom_histogram(binwidth = 2)+
  theme_classic() + 
  theme(axis.text=element_text(size=10))+
  xlab ("Number of iSNV per sample") + ylab ("Number of samples")
  

#isnv frequency (1% bins)
snv_freq_plot <- ggplot (snv_meta, aes (avg_freq))+ geom_histogram(binwidth =0.01)+
  theme_classic ()+ 
  theme(axis.text=element_text(size=10))+
  xlab ("iSNV frequency") + ylab ("Number of iSNV")+
  xlim (0, 0.5)


# divergence rate 
sample_divergence <- read.csv ("Results/SampleDivPerSiteGenomeAllSamples.csv")


mutation.labs <- c("Nonsynonymous", "Synonymous")
names(mutation.labs) <- c("Non", "Syn")


all_sample_plot <- ggplot (DataWithinDiv, aes (x= DPI,y= DivPerDay)) + 
  geom_point (alpha =0.6, shape = 16, aes (color=mutation_type), show.legend = FALSE) +
  scale_color_manual(values=c("#77427F", "#21A386"))+
  geom_crossbar (data = DataWithinDiv_Sum, aes (x = DPI, y =mean, ymin =mean, ymax=mean), linewidth = .4)+
  theme_bw ()+ 
  theme(axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab ("Days post infection") + ylab ("Divergence rate")+
  facet_wrap (vars(mutation_type), labeller = labeller(mutation_type = mutation.labs)) 

## main text figure conjoined 

plot_top <- plot_grid(snv_per_sample_plot, snv_freq_plot,   nrow=1, labels = c ("A", "B"), vjust = 0.5)
plot_bottom <- plot_grid (all_sample_plot, vjust=0.5, labels="C")
plot_grid  ("",plot_top, "",  plot_bottom, nrow = 4, rel_heights = c (.1, 1, .1, 1))




### number of isnv per sample supplemental
isnv_antiviral_plot <- ggplot (data=subset( snv_meta_count), aes (n))+ 
  geom_density( aes (fill = ever_antiviral), alpha =0.5) + 
  theme_bw () + theme(legend.position = "none") +
  ylab ("Density") + xlab ("Number of iSNV per Sample")

isnv_vax_plot <- ggplot (data=subset(snv_meta_count, !is.na(cdc_flu_vx)), aes (n))+ 
  geom_density( aes (fill = cdc_flu_vx), alpha =0.5) + 
  theme_bw () + theme(legend.position = "none") +
  ylab ("Density") + xlab ("Number of iSNV per Sample")

isnv_sex_plot <- ggplot (snv_meta_count, aes (n))+ 
  geom_density( aes (fill = sex), alpha =0.5) + 
  theme_bw () + theme(legend.position = "none") +
  ylab ("Density") + xlab ("Number of iSNV per Sample")


all_isnv_age_plot <- ggplot (snv_meta_count, aes (age.new, n)) + 
  geom_point (size =.3, aes(color = age_class))+
  scale_color_manual("grey", "black")
  theme_cowplot () + xlab ("")+ ylab ("") + 
  theme (panel.background = element_rect(fill = "white"), axis.text = element_text(size=8))


low_isnv_age_plot <- ggplot (snv_meta_count, aes (age.new, n)) + geom_point (size =1)+
  theme_bw () + xlab ("Age") + ylab ("Number of iSNV per Sample") + ylim (0,50)


isnv_time_plot <- ggplot (data=subset( number_of_mutations_per_sample_covariates, !is.na(Days_post_symp)), aes (as.factor (Days_post_symp), n)) + geom_boxplot() +
  theme_cowplot (12) + 
  xlab ("Days Post Symptom Onset") +  ylab ("Number of iSNV") 



plot_top <-  plot_grid ( isnv_vax_plot, isnv_antiviral_plot,
                         labels = c("A", "B"), vjust = 0.5
)
plot_middle  <- plot_grid (isnv_sex_plot, age_plot_inset,
                           align = "h",
                           axis = "bt",
                           labels = c( "C", "D"), vjust = 0.5
)

plot_bottom <- plot_grid (isnv_time_plot_noHigh, labels="E", vjust = 0.5)

plot_grid( NULL,
           plot_top, NULL,
           plot_middle, NULL, plot_bottom,
           nrow=6, rel_heights = c(0.05, 1, 0.05,1, 0.05, .85)
)


### isnv frequency supplemental data 

snv_meta$age_class <- ordered(snv_meta$age_class, levels = c("Child", "Adolescent", "Adult"))
age_freq_plot <-   ggplot (snv_meta, aes (avg_freq)) + geom_density(fill= "grey" ) + 
  facet_grid (rows = vars (age_class)) +
  theme_cowplot (12) +
  scale_x_continuous(expand = c(0, 0), lim = c(0,0.5)) + scale_y_continuous(expand = c(0, 0)) +
  ylab ("Density") + xlab ("iSNV Frequency")

antiviral_freq_plot  <- ggplot (subset( snv_meta, !is.na(ever_antiviral)), aes (avg_freq))+ geom_density(aes (group= ever_antiviral, fill= as.factor (ever_antiviral)), alpha =0.5)+ 
  theme_cowplot (12) + xlim (0,0.5) + 
  ylab ("Density") + xlab ("iSNV Frequency") +
  theme(legend.position = "none") 


vax_freq_plot <-  ggplot (subset( snv_meta, !is.na(cdc_flu_vx)), aes (avg_freq))+ geom_density(aes (group= cdc_flu_vx, fill= cdc_flu_vx), alpha =0.5, kernel ="epanechnikov")+ 
  theme_cowplot (12) + xlim (0,0.5) + 
  ylab ("Density") + xlab ("iSNV Frequency")+ 
  theme(legend.position = "none") 


mutation_type_freq_plot <-  ggplot (snv_meta, aes (avg_freq))+ geom_density(aes (group= mutation_type, fill= mutation_type), alpha =0.5)+ 
  theme_cowplot (12) + xlim (0,0.5) + 
  ylab ("Density") + xlab ("iSNV Frequency") + 
  theme(legend.position = "none") 

days_post_sym_freq_plot <- ggplot (snv_meta, aes ( Days_post_symp, avg_freq)) + geom_jitter (alpha =0.3, shape =16, size =0.7) + 
  ylim (0, 0.5) + theme_cowplot (12)+
  xlab ("Days Post Symptom Onset") + ylab ("iSNV Frequency")

time_freq_sum <- snv_meta %>% group_by (Days_post_symp) %>% summarize (avg_freq = median (avg_freq))

days_post_sym_freq_mini_plot <- 
  ggplot (snv_meta, aes ( Days_post_symp, avg_freq)) + geom_jitter (alpha =0.3, shape =16, size =0.7) + 
  ylim (0, 0.1) + theme_cowplot (12) +
  geom_crossbar(data=time_freq_sum, aes(ymin = avg_freq, ymax = avg_freq),
                col="red", width = .5, size =0.4) +
  xlab ("Days Post Symptom Onset") + ylab ("iSNV Frequency")


days_post_sym_freq_composite_plot <- 
  plot_grid (days_post_sym_freq_plot, NULL,  days_post_sym_freq_mini_plot,  ncol =3,
             labels="E", rel_widths = c (1,0.05,1), vjust = 0.5)

freq_density <- plot_grid (antiviral_freq_plot,  vax_freq_plot,  mutation_type_freq_plot,
                           nrow=3, labels = c ( "B","C", "D"), vjust = 0.5)

age_freq_plot_grid <- plot_grid (NULL, age_freq_plot, nrow =2, labels = c ("A"), vjust = 0.5, rel_heights = c(0.001,1))
top_freq_plot <- plot_grid (age_freq_plot_grid, NULL, freq_density, ncol =3, rel_widths = c (1,0.1,0.8) )

plot_grid (NULL, top_freq_plot, NULL, days_post_sym_freq_composite_plot,
           nrow = 4, ncol = 1, rel_heights =  c (0.05, 1.5, 0.02, .75 ))

