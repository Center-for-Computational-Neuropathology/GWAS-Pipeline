"""
Created on Mon Jan 29 12:22 2024
Goal: Perform SNP and Sample level QC
@@author: Shrishtee Kandoi (github: ShrishteeKandoi)
Username: kandos01
"""

#Load required modules
ml plink
cd /sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/KAT

#Read each chrom per batch
# i=${1}
dir=/sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/KAT

########################################
################ Step 1 ################
########################################

#Investigate missingness per individual and per SNP and make histograms
plink --bfile ${dir}/PART_KatQC_Hg19 --missing

#Generate plots to visualize the missingness results
Rscript --no-save /sc/arion/projects/tauomics/Shrishtee/1_QC_GWAS/hist_miss.R

#Delete SNPs and Samples with missingness > 0.2
plink --bfile ${dir}/PART_KatQC_Hg19 --geno 0.2 --mind 0.2 --make-bed --out ${dir}/PART_KatQC_Hg19.2 -keep-allele-order

#Check missingness again and plot
plink --bfile ${dir}/PART_KatQC_Hg19.2 --missing
Rscript --no-save /sc/arion/projects/tauomics/Shrishtee/1_QC_GWAS/hist_miss.R

########################################
################ Step 2 ################
########################################

# Select autosomal SNPs only (i.e., from chromosomes 1 to 22)
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' ${dir}/PART_KatQC_Hg19.2.bim > snp_1_22.txt
plink --bfile ${dir}/PART_KatQC_Hg19.2 --extract snp_1_22.txt --make-bed --out ${dir}/PART_KatQC_Hg19.3

# Generate a plot of the MAF distribution
plink --bfile ${dir}/PART_KatQC_Hg19.3 --freq --out MAF_check
Rscript --no-save /sc/arion/projects/tauomics/Shrishtee/1_QC_GWAS/MAF_check.R

# Remove SNPs with a low MAF frequency (A conventional MAF threshold for a regular GWAS is between 0.01 or 0.05, depending on sample size)
plink --bfile ${dir}/PART_KatQC_Hg19.3 --maf 0.01 --make-bed --out ${dir}/PART_KatQC_Hg19.4

#MAF - 0.01 ===> ? remained, ? SNPs removed
#MAF - 0.05 ===> ? remained, ? SNPs removed

########################################
################ Step 4 ################
########################################

# Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE) & Check the distribution of HWE p-values of all SNPs
plink --bfile ${dir}/PART_KatQC_Hg19.4 --hardy

awk '{ if ($9 <0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe
Rscript --no-save /sc/arion/projects/tauomics/Shrishtee/1_QC_GWAS/hwe.R
#Error in read.table(file = "plinkzoomhwe.hwe", header = TRUE) : no lines available in input : Execution halted

plink --bfile ${dir}/PART_KatQC_Hg19.4 --hwe 1e-6 --make-bed --out ${dir}/PART_KatQC_Hg19.5
plink --bfile ${dir}/PART_KatQC_Hg19.5 --hwe 1e-10 --hwe-all --make-bed --out ${dir}/PART_KatQC_Hg19.6

# SNPs pass filters and QC

########################################
################ Step 5 ################
########################################

# Inversions file is the region file with high-LD: https://github.com/cran/plinkQC/blob/master/inst/extdata/high-LD-regions-hg19-GRCh37.txt
plink --bfile ${dir}/PART_KatQC_Hg19.6 --extract ../inversions_grch37.txt --range --indep-pairwise 50 5 0.2 --out indepSNP

plink --bfile PART_KatQC_Hg19.6 --extract indepSNP.prune.in --het --out R_check

Rscript --no-save /sc/arion/projects/tauomics/Shrishtee/1_QC_GWAS/check_heterozygosity_rate.R

Rscript --no-save /sc/arion/projects/tauomics/Shrishtee/1_QC_GWAS/heterozygosity_outliers_list.R

sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

plink --bfile PART_KatQC_Hg19.6 --remove het_fail_ind.txt --make-bed --out PART_KatQC_Hg19.7

########################################
################ Step 6 ################
########################################

# Check for relationships between individuals with a pihat > 0.2.
plink --bfile PART_KatQC_Hg19.7 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2

awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome > zoom_pihat.genome


##DONEE



##### STEPS THAT HAVEN'T BEEN PERFORMED


########################################
################ Step 2 ################
########################################

# Check for sex discrepancy.
plink --bfile PART_KatQC_Hg19.7 --check-sex --allow-no-sex

# Generate plots to visualize the sex-check results.
Rscript --no-save gender_check.R

# 1) Delete individuals with sex discrepancy.
grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt
# This command generates a list of individuals with the status �PROBLEM�.
plink --bfile PART_KatQC_Hg19.7 --remove sex_discrepancy.txt --make-bed --out upenn_ucla_mssm.chr$i.6 
# This command removes the list of individuals with the status �PROBLEM�.

# 2) Impute-sex.
#plink --bfile HapMap_3_r3_5 --impute-sex --make-bed --out HapMap_3_r3_6
# This imputes the sex based on the genotype information into your data set.

