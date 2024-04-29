#!/bin/bash

###################################
############### SOW52 ###############
###################################

#Convert plink to VCF
cd /sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/SOW52
plink --bfile SOW52_PART_sex.update --recode vcf --out SOW52_PART_sex_update

#Start a high power node
bsub -P acc_tauomics -q express -n 1 -R rusage[mem=100000] -W 10:00 -Is /bin/bash

#Load necessary libraries
ml crossmap
ml python
ml plink

#Then liftover from hg19 to hg38 using crossmap
# CrossMap.py vcf chain_file hg19_vcf_file hg38_ref_file output_hg38_vcf
cd /sc/arion/projects/tauomics/Shrishtee

CrossMap.py vcf PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz PART_DATA_SOW_KAT_CHOP/SOW52/SOW52_PART_vcf_file.vcf Reference_files/hg38.fa.gz PART_DATA_SOW_KAT_CHOP/SOW52/hg38/SOW52_PART_vcf_hg38.vcf

#Convert vcf file to Plink (BED/BIM/FAM) files
cd /sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/SOW52/hg38
plink --vcf SOW52_PART_vcf_hg38.vcf --make-bed --out SOW52_PART_vcf_hg38 --double-id --allow-extra-chr

###################################
############### KAT ###############
###################################

#Convert plink to VCF

cd /sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/KAT
plink --bfile PART_KatQC_Hg19 --recode vcf --out PART_KatQC_Hg19

#Start a high power node
bsub -P acc_tauomics -q express -n 1 -R rusage[mem=100000] -W 10:00 -Is /bin/bash

#Load necessary libraries

ml crossmap
ml python
ml plink

#Then liftover from hg19 to hg38 using crossmap
# CrossMap.py vcf chain_file hg19_vcf_file hg38_ref_file output_hg38_vcf
cd /sc/arion/projects/tauomics/Shrishtee

CrossMap.py vcf PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz PART_DATA_SOW_KAT_CHOP/KAT/PART_KatQC_Hg19.vcf Reference_files/hg38.fa.gz PART_DATA_SOW_KAT_CHOP/KAT/hg38/PART_KatQC_hg38.vcf

#Convert vcf file to Plink (BED/BIM/FAM) files
cd /sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/KAT/hg38
plink --vcf SOW52_PART_vcf_hg38.vcf --make-bed --out SOW52_PART_vcf_hg38 --double-id --allow-extra-chr

#########LOGS#########

#SOW52
CrossMap.py vcf PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz PART_DATA_SOW_KAT_CHOP/SOW52/SOW52_PART_vcf_file.vcf Reference_files/hg38.fa.gz PART_DATA_SOW_KAT_CHOP/SOW52/SOW52_PART_vcf_hg38
Read the chain file "PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz" 
Filter out variants [reference_allele == alternative_allele] ...
Updating contig field ... 
Lifting over ... 
Total entries: 700078
Failed to map: 106163


#KAT
CrossMap.py vcf PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz PART_DATA_SOW_KAT_CHOP/KAT/PART_KatQC_Hg19.vcf Reference_files/hg38.fa.gz PART_DATA_SOW_KAT_CHOP/KAT/hg38/PART_KatQC_hg38.vcf

Read the chain file "PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz" 
Filter out variants [reference_allele == alternative_allele] ...
Updating contig field ... 
Lifting over ... 
Total entries: 683078
Failed to map: 100661


#Crossmap imputed file
#Michigan server
CrossMap.py vcf PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz PART_DATA_SOW_KAT_CHOP/chr_5/chr5.dose.vcf.gz Reference_files/hg38.fa.gz PART_DATA_SOW_KAT_CHOP/chr_5/chr5.dose_hg38.vcf.gz
CrossMap.py vcf PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz PART_DATA_SOW_KAT_CHOP/chr_6/chr6.dose.vcf.gz Reference_files/hg38.fa.gz PART_DATA_SOW_KAT_CHOP/chr_6/chr6.dose_hg38.vcf.gz

#Minimac
CrossMap.py vcf PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz PART_DATA_SOW_KAT_CHOP/SOW52/Imputed/imputed_PART_SOW52_PART.separate.5.vcf.gz Reference_files/hg38.fa.gz PART_DATA_SOW_KAT_CHOP/SOW52/Imputed/imputed_PART_SOW52.separate.5_hg38.vcf.gz
CrossMap.py vcf PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz PART_DATA_SOW_KAT_CHOP/SOW52/Imputed/imputed_PART_SOW52_PART.separate.6.vcf.gz Reference_files/hg38.fa.gz PART_DATA_SOW_KAT_CHOP/SOW52/Imputed/imputed_PART_SOW52.separate.6_hg38.vcf.gz

###### INSTALL LATEST VERSION OF CROSSMAP ########
mamba install CrossMap
CrossMap vcf PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz PART_DATA_SOW_KAT_CHOP/SOW52/Imputed/imputed_PART_SOW52_PART.separate.6.vcf.gz Reference_files/hg38.fa.gz PART_DATA_SOW_KAT_CHOP/SOW52/Imputed/hg38/imputed_PART_SOW52.separate.6_hg38.vcf.gz


for chr in {1..22}; do
    # Run CrossMap for each chromosome
    CrossMap vcf PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz \
        PART_DATA_SOW_KAT_CHOP/SOW52/Imputed/imputed_PART_SOW52_PART.separate.${chr}.vcf.gz \
        Reference_files/hg38.fa.gz \
        PART_DATA_SOW_KAT_CHOP/SOW52/Imputed/hg38/PART_SOW52_Imputed.chr${chr}_hg38.vcf
    
    # Compress the generated file using bgzip
    bgzip PART_DATA_SOW_KAT_CHOP/SOW52/Imputed/hg38/PART_SOW52_Imputed.chr${chr}_hg38.vcf
    tabix -p vcf PART_DATA_SOW_KAT_CHOP/SOW52/Imputed/hg38/PART_SOW52_Imputed.chr${chr}_hg38.vcf.gz
done

for chr in {1..22}; do
tabix -p vcf PART_SOW52_Imputed.chr${chr}_hg38.vcf.gz
done


2024-04-26 12:13:39 [INFO]  Read the chain file "PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz" 
2024-04-26 12:13:40 [INFO]  Filter out variants [reference_allele == alternative_allele] ...
2024-04-26 12:13:40 [INFO]  Updating contig field ... 
2024-04-26 12:13:40 [INFO]  Lifting over ... 
2024-04-26 12:15:33 [INFO]  Total entries: 2954410
2024-04-26 12:15:33 [INFO]  Failed to map: 1760

# 2024-04-25 05:35:11 [INFO]  Read the chain file "PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz" 
# 2024-04-25 05:35:11 [INFO]  Filter out variants [reference_allele == alternative_allele] ...
# 2024-04-25 05:35:11 [INFO]  Updating contig field ... 
# 2024-04-25 05:35:11 [INFO]  Lifting over ... 
# Traceback (most recent call last):
#   File "/hpc/packages/minerva-centos7/py_packages/3.7/bin/CrossMap.py", line 319, in <module>
#     crossmap_vcf_file(mapping = mapTree, infile= in_file, outfile = out_file, liftoverfile = chain_file, refgenome = genome_file, noCompAllele = args.no_comp_alleles, compress = args.compression, cstyle = args.cstyle)
#   File "/hpc/packages/minerva-centos7/py_packages/3.7/lib/python3.7/site-packages/cmmodule/mapvcf.py", line 183, in crossmap_vcf_file
#     fields[4] = revcomp_DNA(fields[4], True)
#   File "/hpc/packages/minerva-centos7/py_packages/3.7/lib/python3.7/site-packages/cmmodule/utils.py", line 104, in revcomp_DNA
#     return ''.join([complement[base] for base in reversed(seq)])
#   File "/hpc/packages/minerva-centos7/py_packages/3.7/lib/python3.7/site-packages/cmmodule/utils.py", line 104, in <listcomp>
#     return ''.join([complement[base] for base in reversed(seq)])
# KeyError: '>'

#Read the chr6 file to understand by the error is there
library(tidyr)
library(dplyr)
library(data.table)

chr6 <- fread("imputed_PART_SOW52_PART.separate.6.vcf.gz")
subset <- chr6 %>% select ("#CHROM", "POS", "REF", "ALT")

# Check if any column contains '>'
contains_gt <- apply(subset, 2, function(x) any(grepl(">", x)))

# Print columns containing '>'
cols_with_gt <- names(subset)[contains_gt]
if (length(cols_with_gt) > 0) {
  print(cols_with_gt)
} else {
  print("No columns with '>' character.")
}

#ANSWER
# [1] "ALT"

# Print rows that contain '>'
rows_with_gt <- apply(subset, 1, function(row) any(grepl(">", row)))

if (any(rows_with_gt)) {
  print(subset[rows_with_gt, ])
} else {
  print("No rows contain '>' character.")
}

with_symbol <- subset[grep(">", subset$ALT), "ALT"]
unique_alt_with_gt <- unique(subset[grep(">", subset$ALT), "ALT"])

#ANSWER

#      #CHROM       POS REF   ALT
#    1:      6    181655   T <CN2>
#    2:      6    192162   C <CN0>
#    3:      6    209927   G <CN2>
#    4:      6    237754   T <CN0>
#    5:      6    256502   G <CN2>
#   ---                           
# 2424:      6 170533724   T <CN0>
# 2425:      6 170564526   G <CN0>
# 2426:      6 170652358   T <CN0>
# 2427:      6 170760038   T <CN0>
# 2428:      6 170767987   T <CN0>



#VCF-liftover
./VCF-liftover/vcf-liftover.sh PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz PART_DATA_SOW_KAT_CHOP/chr_5/chr5.dose.vcf.gz PART_DATA_SOW_KAT_CHOP/chr_5/chr5.dose_hg38_vcf_liftover.vcf.gz
