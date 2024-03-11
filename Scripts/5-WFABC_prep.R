### prepare files for WFABC 

library (dplyr)
library (tidyr)
library (gdata)
library (stringr)
library(data.table)

## read in SNV file that has the frequency of the ref allele in the first sample for each individual

minor_snv <- read.csv ("Results/SNV_with_corrected_frequency.csv", colClasses = c(household = "character", individual= "character", specimen_number= "character"))
minor_snv <- filter (minor_snv, ALT_FREQ_1 >= 0.005, ALT_FREQ_2 >= 0.005)

## if more than one sample per day, keep self collected sample

minor_wfabc  <-  minor_snv %>%  group_by (hhsubid_pos, Days_post_symp_hh) %>% slice (which.max (collection_type)) 


## keep snv from individuals with 2 or more samples
minor_count <-  minor_wfabc %>% group_by (hhsubid) %>% distinct ( sample, .keep_all=T) %>% count ( hhsubid) %>% filter (n >1)
minor_wfabc  <- filter (minor_wfabc, hhsubid %in% minor_count$hhsubid) 

# get data of first sample for each individual
minor_wfabc$SPECDT_1 <- as.Date (minor_wfabc$SPECDT_1, "%Y-%m-%d")
date_first_sample <- minor_wfabc %>% group_by (hhsubid) %>% slice (which.min (SPECDT_1)) %>%
  select (hhsubid, SPECDT_1) %>% rename (date_first_sample =SPECDT_1)

# separate time of first sample

minor_wfabc <- left_join (minor_wfabc, date_first_sample , by = "hhsubid" ) ## add in date of first sample
minor_wfabc$SPECDT_1 <- as.Date (minor_wfabc$SPECDT_1, "%Y-%m-%d" ) 
minor_wfabc$date_first_sample <- as.Date (minor_wfabc$date_first_sample , "%Y-%m-%d" )
minor_wfabc <- mutate (minor_wfabc, days_since_first_sample = SPECDT_1-date_first_sample ) # calculated days since first sample
minor_wfabc$days_since_first_sample <- as.numeric (minor_wfabc$days_since_first_sample)
minor_wfabc <- mutate (minor_wfabc, gen_since_first_sample = 4* days_since_first_sample) # calculated generations since first sample (assume 6hr generation time)


### add in allele count for each mutation - depth * allele frequency
minor_wfabc <- mutate (minor_wfabc , Total_DP  = TOTAL_DP_1+ TOTAL_DP_2, Allele_count = round (avg_freq_initial_ref* Total_DP))

### seperate by subtype and make input files for each subtype
minor_wfabc$PCR_RESULT_1 <- as.factor (minor_wfabc$PCR_RESULT_1)
for (subtype in levels (minor_wfabc$PCR_RESULT_1)) {
  minor_wfabc_subtype <- filter (minor_wfabc, PCR_RESULT_1 == subtype)
  
  ### get frequency for each timepoint (long form)###

   # allele count df
  minor_wfabc_allele_count <- minor_wfabc_subtype %>% ungroup ()  %>% select ( hhsubid_pos, REGION, POS,  mutation_type, hhsubid,PCR_RESULT_1, onsetdt, age.new, sex, p_seasvx, Year,  hhid,reference, hh_onsetdt, p_seasvx,  initial_mutation, Allele_count, gen_since_first_sample)
  minor_wfabc_allele_count <- minor_wfabc_allele_count %>% group_by (hhsubid_pos) %>% spread (gen_since_first_sample, Allele_count)


  ## separate individuals/mutations with different timepoints

  minor_wfabc_allele_count <-select (minor_wfabc_allele_count,-REGION, -POS, -mutation_type, -hhsubid,-PCR_RESULT_1, -onsetdt, -age.new, -sex, -p_seasvx, -Year, -hhid, -reference, -hh_onsetdt, -initial_mutation)
  minor_wfabc_ac_no_names <-minor_wfabc_allele_count %>% ungroup ( ) %>%  select (-hhsubid_pos)

 
# get combinations of time points  
  make_combinations <- function(x) {

    l <- length(x)
    mylist <- lapply(2:l, function(y) {
      combn(x, y, simplify = FALSE)
    })
    mylist
    unlist(mylist, recursive = FALSE)
  }
}
  results <- make_combinations(colnames(minor_wfabc_names))

  ### Make individual input files for the different time point combinations

  for (i in 1:length(results)) {
    NA_names <- minor_wfabc_allele_count %>% ungroup %>% select ( -results [[i]], -hhsubid_pos)
    NA_names <- list(colnames(NA_names)) %>% unlist (.)
    
    
    if (length (NA_names)==0) {
      ac.df <- minor_wfabc_allele_count %>% filter_at(c(results[[i]]), all_vars(!is.na(.))) %>% 
        select (hhsubid_pos, results[[i]]) 
    }else {
      ac.df <- minor_wfabc_allele_count %>% filter_at(c(results[[i]]), all_vars(!is.na(.))) %>% 
        filter_at(c(NA_names), all_vars(is.na(.)))%>% 
        select (hhsubid_pos, results[[i]]) %>% ungroup ()
        }
    
    if (dim(ac.df)[1] != 0) {
      wfabc_ac.df <- ac.df %>% select(2:ncol(.) ) # get frequency data frame
      depth.df <- dp.df %>% select(2:ncol(.) )  # get coverage data frame
      names <-  ac.df %>% select (1) # get names of individuals and use to make csv file
      generation_col <- colnames (wfabc_ac.df)  # get column names to use in  output files
      generation_name <- paste(generation_col,collapse="_")
    
  

      file_name <- paste ("Results/WFABC/", subtype, "/Actual_depth/", generation_name, "_names.csv", sep="")
      write.csv(names, file_name, row.names = F, quote=F)
    
      header_info <- c (nrow (wfabc_ac.df), ncol (wfabc_ac.df))# get header info to add to top of csv file
      header_file <- paste ("Results/WFABC/", subtype, "/Actual_depth/", generation_name, "_header.txt", sep="")
      write(header_info, header_file)
    
      depth.df <- as.matrix(depth.df)
      wfabc_ac.df <- as.matrix(wfabc_ac.df)
      wfabc_final <- interleave (depth.df, wfabc_ac.df) 
      wfabc_file <- paste ("Results/WFABC/",subtype, "/Actual_depth/", generation_name, ".txt", sep="")
      write.csv (wfabc_final, wfabc_file, row.names=F, quote=F)
    }
  }
  
  

  # make input file for Ne estimation- use first 2 columns if time between generations is 4
  # make empty data frame (one for depth )
  Ne_ac.df <-data.frame(Col1 = character(),
                   Col2 = character())
  Ne_ac.df <- rename (Ne_ac.df, "0" = 1, "4" =2)

 
  # make input file using allele counts and actual depth 
  for (row in 1:nrow(minor_wfabc_ac_no_names)){
    row_wfabc <- minor_wfabc_ac_no_names [row, ]
    row_wfabc <-row_wfabc %>% select_if(~ !any(is.na(.)))
    
  
    if (ncol (row_wfabc) > 1) {
      my_colnames <- as.numeric (colnames(row_wfabc))
      if (my_colnames[2]-my_colnames [1] ==4){
        row_wfabc <- row_wfabc %>% select ( 1,2) %>% rename ("0" = 1, "4" =2)
        row_wfabc [,1]<- as.character(row_wfabc [,1])
        row_wfabc [,2]<- as.character(row_wfabc [,2])
        Ne_ac.df <- bind_rows ( Ne_ac.df ,row_wfabc )
      }
    }
  }
  

  Ne_final_ac <- interleave (Ne_dp.df, Ne_ac.df )
  Ne_file <- paste ("Results/WFABC/", subtype, "/Actual_depth/Ne_sample_size_input.txt", sep="")
  write.csv (Ne_final_ac , Ne_file , row.names=F, quote=F)

  
  Ne_af.df <-data.frame(Col1 = character(),
                        Col2 = character())
  Ne_af.df <- rename (Ne_af.df, "0" = 1, "4" =2)
  



  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  