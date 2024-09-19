FluTES Within Host Evolution
--
This repository contains code and small intermediate data files associated with the manuscript "Influenza A virus within-host evolution and positive selection in a densely sampled household cohort over three seasons". 

The directory is organized as follows:

**Pipeline** -  This folder contains shell scripts for processing raw fastq files and calling variants. See manuscript for details of sample processing .  Raw sequencing data is available in the SRA as BioProject PRJNA1085292.  

**Processed_data/Secondary_processing** This file contains the variant file generated from ivar with associated metadata. The variants have been fully filtered and is the input file for all downstram analyses

**References** - This folder contains the reference sequence and gff file for each seasons (2017/18 - 2019/20) vaccine strain. Note that the references include the universal IAV priming sites. Downstream variant files use these reference coordinants. It also includes files used in the divergence rate calculations.  

**Metadata** - This folder contains metadata for all specimens. The "Definitions.xlsx" file contains the definitions for the variable columns. 

**Scripts** - This folder contains custom R scripts to perform analyses in the manuscript. The files are numbered in order of usage with figure scripts being run last. Most intermediate files are included in this repo and each script should be able to run independently of each other. 

**Results** - This folder contains intermediate files and results of analyses. 
