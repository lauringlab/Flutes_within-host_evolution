
## Allele Trajectoyry- looking for an increase in frequency

#Keep mutations that are in 2 or more samples withion an individual
  # if 2 samples were collected on the same day keep the self-collected sample

# load libraries 
library (dplyr)
library (tidyr)
library (ggplot2)

# read in data
IAV <- read.csv ("Processed_data/Secondary_processing/SNV_with_meta_data.csv", 
                 colClasses = c(household = "character", individual= "character", specimen_number= "character"))
IAV$SPECDT_1 <- as.Date(IAV$SPECDT_1, "%m/%d/%y")

## get SNV that are in more than 1 sample within a person
 Position_count <- count (IAV, REGION, POS, hhsubid ) %>%
   filter ( n>1) %>%
   unite ( hhsubid_pos, c (hhsubid,REGION, POS), sep= "-", remove =F)
 
 mult_position <- IAV %>% unite ( hhsubid_pos, c (hhsubid,REGION, POS), sep= "-", remove =F) %>%
  filter ( hhsubid_pos %in%  Position_count$hhsubid_pos)

# check if any mutations changed from the minor allele to the major allele
mult_mutations_count <-  mult_position %>% distinct ( mutation, hhsubid, .keep_all=T) %>%
  count (hhsubid_pos) %>% filter (n>1)

mult_mutations <-  mult_position %>% filter (hhsubid_pos %in% mult_mutations_count$hhsubid_pos)

mult_mutations_count1 <-  mult_position %>% distinct ( mutation, hhsubid, .keep_all=T) %>%
  count (hhsubid_pos) %>% filter (n==1)

single_mut <- mult_position %>% filter (hhsubid_pos %in%  mult_mutations_count1$hhsubid_pos) # data frame with one mutation per position


# check for more than 2 alleles 
plus_2_alleles <- mult_mutations %>% gather (allele_type, allele, REF, ALT) %>% 
  distinct (hhsubid_pos, allele) %>% count (hhsubid_pos) %>% filter (n>2)


  
mult_mutations <-filter (mult_mutations, !(hhsubid_pos %in% plus_2_alleles$hhsubid_pos)) # filter out positions that have 3 or more alleles

# use first sample as the REF allele for all samples with mutations at the same position

mult_mutations$SPECDT_1 <- as.Date(mult_mutations$SPECDT_1, "%m/%d/%y")
first_mutation <- mult_mutations %>% group_by (hhsubid_pos) %>% arrange (mult_mutations$SPECDT_1) %>% 
  summarise (first (mutation)) %>% rename (initial_mutation = 2) # get first mutation

mult_mutations <- full_join (mult_mutations, first_mutation, by = "hhsubid_pos") 

mult_mutations <- mutate (mult_mutations, avg_freq_initial_ref = ifelse (mutation == initial_mutation, avg_freq, 1-avg_freq))  
write.csv (mult_mutations, "Results/Major_minor_allele_trajectory.csv", row.name=F, quote=F)

## add back in positions with a single mutation for wfabc analysis and for allele trajectory of all alleles

single_mut <- mutate (single_mut,  avg_freq_initial_ref = avg_freq, initial_mutation = mutation)
allele_traj.df <- bind_rows( single_mut, mult_mutations)
write.csv (allele_traj.df, "Results/SNV_with_corrected_frequency.csv", row.names=F, quote=F)

## add in positions with a single sample 
Position_count_1 <- count (IAV, REGION, POS, hhsubid ) %>%
  filter ( n == 1) %>%
  unite ( hhsubid_pos, c (hhsubid,REGION, POS), sep= "-", remove =F)

single_position <- IAV %>% unite ( hhsubid_pos, c (hhsubid,REGION, POS), sep= "-", remove =F) %>%
  filter ( hhsubid_pos %in%  Position_count_1$hhsubid_pos)

single_position <- mutate (single_position,  avg_freq_initial_ref = avg_freq, initial_mutation = mutation)
Corrected_freq_all.df <- bind_rows( single_position,allele_traj.df )
write.csv (Corrected_freq_all.df , "Results/SNV_with_corrected_frequency_all_mutations.csv", row.names=F, quote=F)


## iSNV that go from major to minor 

Major_minor_list <- mult_mutations %>% filter (initial_mutation != mutation) # get list of iSNV that go from minor to major
Major_minor <- filter (mult_mutations, hhsubid_pos %in% Major_minor_list$hhsubid_pos)

Major_minor_single_date <-Major_minor %>% distinct (hhsubid_pos, SPECDT_1) %>% 
  count (hhsubid_pos) %>% filter (n >1) # only include individuals with iSNV on multiple days
Major_minor <- filter (Major_minor, hhsubid_pos %in% Major_minor_single_date$hhsubid_pos)

write.csv (Major_minor, "Results/Major_minor_allele_trajectory.csv", row.name=F, quote=F)


# get iSNV that go from minor to major in antigenic sites
Major_minor_H3N2 <- Major_minor %>% filter (PCR_RESULT_1 == "FluAH3", REGION == "HA", POS %in% antigenic_sites$H3N2) 
write.csv (Major_minor_H3N2 , "Results/H3N2_antigenic_sites_major_minor_allele_trajectory.csv", row.names=F, quote=F) 


Major_minor_H1N1 <- Major_minor %>% filter (PCR_RESULT_1 == "FluAH1N1", REGION == "HA", POS %in% antigenic_sites$H1N1) 

