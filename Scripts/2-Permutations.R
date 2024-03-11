library(tidyverse)
library(wesanderson)

set.seed(42)
palette <- wes_palette("Darjeeling1")

snv_meta <- read.csv ("Processed_data/Secondary_processing/SNV_with_meta_data.csv", colClasses = c(household = "character", individual= "character", specimen_number= "character"))





# for iSNV that are in multiple individuals in a household, only keep 1 instance of iSNV
variants <- snv_meta  %>% distinct (hhid, mutation, .keep_all=T)
variants <- filter (variants, Year== c ("17"),  PCR_RESULT_1 == "FluAH3") # do for every combination of Year and subtype (vaccine reference) 
# ========================== How many mutations to draw per individual? ====================
mutationsPerIndividual <- variants  %>% group_by(hhsubid) %>%
  
  summarize(numMutations = length(unique(mutation)))  # distribution to draw from



numIndividuals <- length(unique(variants$hhsubid)) # number of individuals to simulate

# =========================== Support functions ===========================

single_permutation <- function(sites)
  
{
  
  sampled_sites <- c()
  
  for(index in 1:numIndividuals)
    
  {
    
    mutationsPerIndividual_drawn <- sample(mutationsPerIndividual$numMutations, size = 1)
    
    sites_one_individual <- sample(sites, size = mutationsPerIndividual_drawn, replace = FALSE)
    
    sampled_sites <- c(sampled_sites, sites_one_individual)
    
  }
  
  
  
  frequencies <- data.frame(table(sampled_sites)) %>% filter(Freq > 1)
  
  
  
  shared_2 <- nrow(filter(frequencies, Freq == 2))
  
  shared_3 <- nrow(filter(frequencies, Freq == 3))
  
  shared_4 <- nrow(filter(frequencies, Freq == 4))
  
  shared_high <- nrow(filter(frequencies, Freq > 4))
  
  
  
  mutations_shared <- data.frame(shared2 = shared_2, shared3 = shared_3, shared4 = shared_4, sharedHigh = shared_high)
  
  
  
  return(mutations_shared)
  
}


# ======================= Function to run one test =====================


run_permutation_test <- function(permutations = 1000, fraction_mutable)
  
{
  
  
  
  possible_sites <- 13101 #based on H1N1 michigan orf size ( Brisbane- 13125, Hong Kong-13101, Singapore-13107)
  
  results <- data.frame()
  
  for(index in 1:permutations)
    
  {
    
    permutation_result <- single_permutation(possible_sites)
    
    results <- rbind(results, permutation_result)
    
  }
  
  
  
  return(results)
  
}
# ============================ Plot observed mutations vs. null distribution ==========================




results_plot<- run_permutation_test(fraction_mutable = 0.6)

num_shared_2 <- median(results_plot$shared2)

num_shared_3 <- median(results_plot$shared3)

num_shared_4 <- median(results_plot$shared4)

num_shared_high <- median(results_plot$sharedHigh)




mutation_by_individual <- distinct (variants, hhsubid, mutation, mutation_type, REGION, POS)
var_count <- count (mutation_by_individual, mutation, mutation_type, REGION, POS, sort=T)
var_count_multiple <- filter (var_count, n>1)
var_count_multiple <- mutate (var_count_multiple, group = ifelse (n ==2, "2", ifelse (n==3, "3", ifelse(n==4, "4", ">4"))))
var_count_multiple$group <- factor(var_count_multiple$group, levels = c("2", "3", "4", ">4"))





# Plot observed and simulated next to each other on the same axis.

var_count_multiple <- mutate(var_count_multiple, type = "Observed")

results_plot <- mutate(results_plot, type = "Simulated")

shared.mutations.null <- ggplot() + 
  
  geom_bar(data = var_count_multiple, aes(group)) +
  
  ylab("Number of Mutations") + 
  
  xlab("Number of Individuals") + 
  
  theme_bw() +
  
  geom_boxplot(data = results_plot, aes(x = "2", y = shared2), alpha = 1) +
  
  geom_boxplot(data = results_plot, aes(x = "3", y = shared3), alpha = 1) +
  
  geom_boxplot(data = results_plot, aes(x = "4", y = shared4), alpha = 1) +
  
  geom_boxplot(data = results_plot, aes(x = ">4", y = sharedHigh), alpha = 1) +
  
  theme(legend.position = "left") +
  
  facet_wrap(~type)+
  ggtitle ("H3N2 2017")

# Observed numbers of shared mutations

var_count_multiple %>% group_by(group) %>% summarize(num = n()) -> observed_counts


observed_counts_num <- observed_counts$num


num_observed_2 <- observed_counts_num[1]

num_observed_3 <- observed_counts_num[2]

num_observed_4 <- observed_counts_num[3]


num_observed_high <- observed_counts_num[4]


# ========================= Get p-values at each fraction =====================

# p-value is the fraction of permutations with a number of shared mutations at or above what is observed for that group.



fractions_to_test <- 0.6 #seq(0.1, 1, by = 0.1)

pvalues_all <- data.frame()

for(fraction in fractions_to_test)
  
{
  
  print(paste0("On fraction ", as.character(fraction)))
  
  test_result <- run_permutation_test(fraction_mutable = fraction)
  
  
  
  pvalues <- c(nrow(filter(test_result, shared2 >= num_observed_2)) / nrow(test_result), 
               
               nrow(filter(test_result, shared3 >= num_observed_3)) / nrow(test_result), 
               
               nrow(filter(test_result, shared4 >= num_observed_4)) / nrow(test_result), 
               
               nrow(filter(test_result, sharedHigh >= num_observed_high)) / nrow(test_result) )
  
  groups <- c("2", "3", "4", ">4")
  
  pvalues_singlefraction <- data.frame(frac = fraction, group = groups, pvalue = pvalues)
  
  pvalues_all <- rbind(pvalues_all, pvalues_singlefraction)
  
}





pval_palette <- wes_palette("Zissou1")

pvalues_all$group <- factor(pvalues_all$group, levels = c("2", "3", "4", ">4"))




#### plot all subtypes
top_plot <- plot_grid(
  shared.mutations.null_Michigan, 
  shared.mutations.null_Brisbane,
  labels = c("A", "B"), vjust = 0.25
)

bottom_plot <- plot_grid( shared.mutations.null_HongKong, 
                          shared.mutations.null_Singapore,
                          labels = c("C", "D"), vjust = 0.25)

plot_grid(NULL, top_plot, NULL, bottom_plot, nrow=4, rel_heights = c(.1, 1,.1, 1))
  

