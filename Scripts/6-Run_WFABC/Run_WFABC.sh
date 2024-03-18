

## get Ne estimates
./wfabc_1 Results/WFABC/FluH1N1/Ne_sample_size_input.txt
./wfabc_1 Results/WFABC/FluH3/Ne_sample_size_input.txt


# run WFABC for each subtype and Ne value
Snakemake -s Scripts/6-Run_WFABC/selection_snakemake_H1N1_Mean_Ne.txt --cores 4
Snakemake -s Scripts/6-Run_WFABC/selection_snakemake_H1N1_Pos_95_Ne.txt --cores 4
Snakemake -s Scripts/6-Run_WFABC/selection_snakemake_H1N1_Neg_95_Ne.txt --cores 4

Snakemake -s Scripts/6-Run_WFABC/selection_snakemake_H3N2_Mean_Ne.txt --cores 4
Snakemake -s Scripts/6-Run_WFABC/selection_snakemake_H3N2_Pos_95_Ne.txt --cores 4
Snakemake -s Scripts/6-Run_WFABC/selection_snakemake_H3N2_Neg_95_Ne.txt --cores 4

# compile results
Snakemake -s Scripts/6-Run_WFABC/merge_wfabc_runs_snakemake.txt --cores 1