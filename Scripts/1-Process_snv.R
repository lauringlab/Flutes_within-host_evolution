library (dplyr)
library (stringr)
library (patchwork)

###### Filter SNV #######

###Coverage
# read avg coverage and merge into one file
avg_cov_HongKong <- read.table ("Processed_data/HongKong14_A/AvgCoverage.all", header=T)
avg_cov_Michigan  <- read.table ("Processed_data/Michigan15_A/AvgCoverage.all", header=T)
avg_cov_Brisbane <-  read.table ("Processed_data/Brisbane18_A/AvgCoverage.all", header=T)
avg_cov_Singapore <- read.table ("Processed_data/Singapore16_A/AvgCoverage.all", header=T)

avg_cov <- bind_rows ( avg_cov_HongKong, avg_cov_Michigan, avg_cov_Brisbane, avg_cov_Singapore )

avg_cov <- separate (avg_cov, sample, c ("sample", "replicate"), "_") # split sample column
avg_cov <- spread (avg_cov, replicate, mean)
avg_cov<- rename (avg_cov, rep_2= 3, rep_1= 2)
avg_cov <- mutate (avg_cov, Pass_coverage = ifelse (rep_1 >= 1000 & rep_2 >=1000, "Pass", "Fail")) # filter by coverage 

write.csv (avg_cov, "Processed_data/Secondary_processing/Pass_coverage.csv", quote=F, row.names=F)

avg_cov_pass <- filter (avg_cov, Pass_coverage == "Pass")

# read in snv 
SNV_HongKong <- read.table ("Processed_data/HongKong14_A/all_variants_filtered", header=T, stringsAsFactors=F)
SNV_Michigan <- read.table ("Processed_data/Michigan15_A/all_variants_filtered", header=T, stringsAsFactors=F)
SNV_Brisbane <- read.table ("Processed_data/Brisbane18_A/all_variants_filtered", header=T, stringsAsFactors=F)
SNV_Singapore <- read.table ("Processed_data/Singapore16_A/all_variants_filtered", header=T, stringsAsFactors=F)


Brisbane_length <- read.csv ("metadata/Brisbane_length_to_add.csv")
  Brisbane_length <- rename (Brisbane_length, REGION=Seg)
Singapore_length <- read.csv ("metadata/Singapore_length_to_add.csv")
  Singapore_length <- rename (Singapore_length, REGION=Seg)
HongKong_length <- read.csv ("metadata/HongKong_length_to_add.csv")
  HongKong_length <- rename (HongKong_length, REGION=Seg)
Michigan_length <- read.csv ("metadata/Michigan_length_to_add.csv")
  Michigan_length <- rename (Michigan_length, REGION=Seg)

SNV_HongKong <- SNV_HongKong %>% left_join (HongKong_length, by== "REGION" ) %>% mutate (concat.pos = POS + Length_to_add)
SNV_Brisbane <- SNV_Brisbane%>% left_join (Brisbane_length, by== "REGION" ) %>% mutate (concat.pos = POS + Length_to_add)
SNV_Michigan <- SNV_Michigan %>% left_join (Michigan_length, by== "REGION" ) %>% mutate (concat.pos = POS + Length_to_add)
SNV_Singapore<- SNV_Singapore%>% left_join (Singapore_length, by== "REGION" ) %>% mutate (concat.pos = POS + Length_to_add)

snv_meta <- bind_rows (SNV_Brisbane, SNV_HongKong)
snv_meta <- bind_rows (snv_meta, SNV_Michigan)
snv_meta <- bind_rows (snv_meta, SNV_Singapore)


# combine and filter samples by GFF 

SNV <- bind_rows (SNV_HongKong, SNV_Michigan, SNV_Brisbane,SNV_Singapore )
SNV<- filter (SNV, sample %in% avg_cov_pass$sample )

SNV_filter_GFF <- filter (SNV, str_detect (GFF_FEATURE, REGION)) # need this later for synonynous/nonsynymous, might only keep coding regions (?)



####----Mutation Type -----###

SNV_filter_GFF <- unite (SNV_filter_GFF, sample_mutation, c (sample, mutation), sep= "-", remove =F)
SNV_filter_GFF<- mutate (SNV_filter_GFF , mutation_type=ifelse(REF_AA==ALT_AA, "Syn", "Non"))

#get mutations in multiple ORF
Count_mutations <- count (SNV_filter_GFF, sample_mutation)
multiple_orf_list <- filter (Count_mutations, n >1)
multiple_orf <- filter (SNV_filter_GFF, sample_mutation %in% multiple_orf_list$sample_mutation )
multiple_orf_select <- select (multiple_orf, GFF_FEATURE, REF_AA, ALT_AA, sample_mutation, mutation, sample, mutation_type)

# change to one mutation per line and get synonymous or nonsynonymous
multiple_orf_distinct <- distinct (multiple_orf_select, sample_mutation, .keep_all=T)
multiple_orf_distinct2<- anti_join (multiple_orf_select, multiple_orf_distinct, by = c("sample_mutation", "GFF_FEATURE"))
multiple_orf_wide <- full_join (multiple_orf_distinct,multiple_orf_distinct2, by = c ("sample_mutation" , "sample", "mutation"))
multiple_orf_wide <- mutate (multiple_orf_wide , mutation_type = ifelse (mutation_type.x =="Non" | mutation_type.y== "Non", "Non", "Syn"))
multiple_orf_mutation_type <- select(multiple_orf_wide, sample_mutation, mutation_type)

multiple_orf <- select (multiple_orf, -mutation_type)
multiple_orf <- full_join (multiple_orf, multiple_orf_mutation_type, by = "sample_mutation")
multiple_orf_distinct_mutation <- distinct (multiple_orf, sample_mutation, .keep_all=T)

#mutations in single ORF
single_orf <- filter (SNV_filter_GFF, ! (sample_mutation %in% multiple_orf_list$sample_mutation ))
all_orf_mutation_type <- bind_rows (single_orf , multiple_orf_distinct_mutation )

# add in average frequency
all_orf_mutation_type <- mutate (all_orf_mutation_type, avg_freq = (ALT_FREQ_1+ ALT_FREQ_2)/2)
write.csv (all_orf_mutation_type, "Processed_data/Secondary_processing/SNV_with_mutation_type.csv", quote=F, row.names=F)

## add in meta_data
## read in data and add meta data 
snv <- read.csv ("Processed_data/Secondary_processing/SNV_with_mutation_type.csv")
meta <- read.csv ("metadata/Vanderbilt_metadata_all_years.csv", colClasses = c(household = "character", individual= "character", specimen_number= "character"))

meta <- rename (meta, "sample"= "CaseID")
snv_meta <- left_join (snv, meta, by ="sample")
snv_meta <- mutate (snv_meta, avg_freq= (ALT_FREQ_1+ALT_FREQ_2)/2)
write.csv (snv_meta , "Processed_data/Secondary_processing/SNV_with_meta_data.csv", row.names=F, quote=F )

##check high number of iSNV
high_snv_samples <- filter  (snv_meta_count, n >= 50  )
high_snv <- filter (snv_meta_count, Study_ID %in% high_snv_samples$Study_ID)
write.csv (high_snv_samples,"Results/high_snv_samples.csv", row.names=F, quote =F)








