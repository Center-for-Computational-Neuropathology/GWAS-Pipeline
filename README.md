# GWAS_Pipeline_CraryLab

Date: April 22 2024

Principal Investigator: Kurt Farrell
Authors: Shrishtee Kandoi

## STEPS INCLUDED IN THE PIPELINE

1. Sample and SNP QC: Data Cleaning
2. Imputation
3. Relatedness
4. Population stratification
5. Run GWAS
6. Visualize results on Locus Zoom
7. Run case-case GWAS
8. Calculate Polygenic Risk Scores / Pathway scores


### 1. QC (Sample and SNP) - Data Clean

```

#Load required modules
ml plink
cd /sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/KAT

#Read each chrom per batch
# i=${1}
dir=/sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/KAT
```
########################################
################ Step 1 ################
########################################

```
#Investigate missingness per individual and per SNP and make histograms
plink --bfile ${dir}/PART_KatQC_Hg19 --missing

#Generate plots to visualize the missingness results
Rscript --no-save /sc/arion/projects/tauomics/Shrishtee/1_QC_GWAS/hist_miss.R

#Delete SNPs and Samples with missingness > 0.2
plink --bfile ${dir}/PART_KatQC_Hg19 --geno 0.2 --mind 0.2 --make-bed --out ${dir}/PART_KatQC_Hg19.2 -keep-allele-order

#Check missingness again and plot
plink --bfile ${dir}/PART_KatQC_Hg19.2 --missing
Rscript --no-save /sc/arion/projects/tauomics/Shrishtee/1_QC_GWAS/hist_miss.R
```

########################################
################ Step 2 ################
########################################

```

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

```

########################################
################ Step 4 ################
########################################

```
# Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE) & Check the distribution of HWE p-values of all SNPs
plink --bfile ${dir}/PART_KatQC_Hg19.4 --hardy

awk '{ if ($9 <0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe
Rscript --no-save /sc/arion/projects/tauomics/Shrishtee/1_QC_GWAS/hwe.R
#Error in read.table(file = "plinkzoomhwe.hwe", header = TRUE) : no lines available in input : Execution halted

plink --bfile ${dir}/PART_KatQC_Hg19.4 --hwe 1e-6 --make-bed --out ${dir}/PART_KatQC_Hg19.5
plink --bfile ${dir}/PART_KatQC_Hg19.5 --hwe 1e-10 --hwe-all --make-bed --out ${dir}/PART_KatQC_Hg19.6

# SNPs pass filters and QC

```

########################################
################ Step 5 ################
########################################

```

# Inversions file is the region file with high-LD: https://github.com/cran/plinkQC/blob/master/inst/extdata/high-LD-regions-hg19-GRCh37.txt
plink --bfile ${dir}/PART_KatQC_Hg19.6 --extract ../inversions_grch37.txt --range --indep-pairwise 50 5 0.2 --out indepSNP

plink --bfile PART_KatQC_Hg19.6 --extract indepSNP.prune.in --het --out R_check

Rscript --no-save /sc/arion/projects/tauomics/Shrishtee/1_QC_GWAS/check_heterozygosity_rate.R

Rscript --no-save /sc/arion/projects/tauomics/Shrishtee/1_QC_GWAS/heterozygosity_outliers_list.R

sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

plink --bfile PART_KatQC_Hg19.6 --remove het_fail_ind.txt --make-bed --out PART_KatQC_Hg19.7

```

########################################
################ Step 6 ################
########################################

```

# Check for relationships between individuals with a pihat > 0.2.
plink --bfile PART_KatQC_Hg19.7 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2

awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome > zoom_pihat.genome


```



##### Optional Step to check for sex discrepancy

```

# Check for sex discrepancy.
plink --bfile PART_KatQC_Hg19.7 --check-sex --allow-no-sex

# Generate plots to visualize the sex-check results.
Rscript --no-save gender_check.R

1. Delete individuals with sex discrepancy.

# This command generates a list of individuals with the status �PROBLEM�.
grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt

# This command removes the list of individuals with the status �PROBLEM�.
plink --bfile PART_KatQC_Hg19.7 --remove sex_discrepancy.txt --make-bed --out upenn_ucla_mssm.chr$i.6 

# 2) Impute-sex.
#plink --bfile HapMap_3_r3_5 --impute-sex --make-bed --out HapMap_3_r3_6
# This imputes the sex based on the genotype information into your data set.

```

## Preparation of data prior to submission to Imputation server


```
#LOAD REQUIRED MODULES

ml plink
ml vcftools
ml bcftools #To use bgzip
ml minimac4

dir=/sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/KAT

#Make the .frq file
plink --bfile ${dir}/PART_KatQC_Hg19.7 --keep-allele-order --freq --out ${dir}/PART_KatQC_Hg19.7 --allow-no-sex

#Separate into chromosome
for chr in {1..22} X Y; do
    # Extract variants for this chromosome
    plink --bfile ${dir}/PART_KatQC_Hg19.7 --chr $chr --make-bed --out ${dir}/chrom/PART_KatQC_Hg19_separate.$chr
done

#Make vcf files
# plink --bfile ${dir}/PART_KatQC_Hg19.7 --recode vcf --out ${dir}/PART_KatQC_Hg19.7 --keep-allele-order

for chnum in {1..22};
  do
  plink --bfile ${dir}/chrom/PART_KatQC_Hg19_separate.$chnum --recode vcf --chr $chnum --out ${dir}/chrom/PART_KatQC_Hg19_separate.$chnum --allow-no-sex
done

#Then sort and zip
# vcf-sort ${dir}/PART_KatQC_Hg19.7.vcf | bgzip -c > pre_impute_PART_KatQC_Hg19.7.vcf.gz

#BGZIP vcf files
for chnum in {1..22};
  do
  vcf-sort ${dir}/chrom/PART_KatQC_Hg19_separate.$chnum.vcf | bgzip -c > ${dir}/chrom/pre_impute_PART_KatQC_Hg19_separate.$chnum.vcf.gz
done

#Index vcf files
for chnum in {1..22};
  do
  tabix -p vcf ${dir}/chrom/pre_impute_PART_KatQC_Hg19_separate.$chnum.vcf.gz
done

#and then you are ready to submit to the imputation server or use minimac4/1.6 (latest version)
# minimac4 1000g_phase3_v5.chr14.with_parameter_estimates.msav pre_impute_upenn_ucla_mssm_impute_chr14.output_file.vcf.gz -o imputed_upenn_ucla_mssm_impute_chr14.output_file.vcf.gz

msav_files=/sc/arion/projects/tauomics/Shrishtee/Imputation/g1k_p3_msav_files_with_estimates

for chnum in {1..22};
  do
  minimac4 ${msav_files}/1000g_phase3_v5.chr$chnum.with_parameter_estimates.msav ${dir}/chrom/pre_impute_PART_KatQC_Hg19_separate.$chnum.vcf.gz -o ${dir}/Imputed/imputed_PART_KatQC_Hg19_separate.$chnum.vcf.gz
done
```


# 2. Imputation

a. Michigan Imputation Server

https://imputationserver.sph.umich.edu/index.html

<!-- Build: hg19
Reference Panel: 1000g-phase-3-v5 (hg19)
Population: mixed
Phasing: eagle
Mode: imputation -->

or 

b. Minimac4: https://github.com/statgen/Minimac4

<!-- Reference panel: https://share.sph.umich.edu/minimac4/panels/ -->

```
ml minimac4

i=${1}

dir=/sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/KAT
msav_files=/sc/arion/projects/tauomics/Shrishtee/Imputation/g1k_p3_msav_files_with_estimates

minimac4 ${msav_files}/1000g_phase3_v5.chr$i.with_parameter_estimates.msav ${dir}/chrom/pre_impute_PART_KatQC_Hg19_separate.$i.vcf.gz -o ${dir}/Imputed/imputed_PART_KatQC_Hg19_separate.$i.vcf.gz

```
