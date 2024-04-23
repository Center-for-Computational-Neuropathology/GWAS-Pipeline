# GWAS_Pipeline_CraryLab

Date: April 22, 2024


* Principal Investigator: Kurt Farrell
* Authors: Shrishtee Kandoi

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

```bash
#Load required modules
ml plink
cd /sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/KAT

#Read each chrom per batch
# i=${1}
dir=/sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/KAT
```

#### a. STEP 1: Missingness of SNPs and individuals

```bash
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


##### Optional Step to check for sex discrepancy

```bash

# Check for sex discrepancy.
plink --bfile PART_KatQC_Hg19.7 --check-sex --allow-no-sex

# Generate plots to visualize the sex-check results.
Rscript --no-save gender_check.R

1. Delete individuals with sex discrepancy.

# This command generates a list of individuals with the status "PROBLEM".
grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt

# This command removes the list of individuals with the status "PROBLEM".
plink --bfile PART_KatQC_Hg19.7 --remove sex_discrepancy.txt --make-bed --out PART_KatQC_Hg19.8

# 2) Impute-sex.
#plink --bfile PART_KatQC_Hg19.8 --impute-sex --make-bed --out PART_KatQC_Hg19.9
# This imputes the sex based on the genotype information into your data set.

```

#### b. STEP 2: Minor allele frequency

```bash

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

#### c. STEP 3: Hardy-Weinberg equilibrium (HWE)

```bash
#Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE) & Check the distribution of HWE p-values of all SNPs
plink --bfile ${dir}/PART_KatQC_Hg19.4 --hardy

awk '{ if ($9 <0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe
Rscript --no-save /sc/arion/projects/tauomics/Shrishtee/1_QC_GWAS/hwe.R

plink --bfile ${dir}/PART_KatQC_Hg19.4 --hwe 1e-6 --make-bed --out ${dir}/PART_KatQC_Hg19.5
plink --bfile ${dir}/PART_KatQC_Hg19.5 --hwe 1e-10 --hwe-all --make-bed --out ${dir}/PART_KatQC_Hg19.6

# SNPs pass filters and QC

```

#### d. STEP 4: Heterozygosity

```bash
# Inversions file is the region file with high-LD: https://github.com/cran/plinkQC/blob/master/inst/extdata/high-LD-regions-hg19-GRCh37.txt

plink --bfile ${dir}/PART_KatQC_Hg19.6 --extract ../inversions_grch37.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
plink --bfile PART_KatQC_Hg19.6 --extract indepSNP.prune.in --het --out R_check
Rscript --no-save /sc/arion/projects/tauomics/Shrishtee/1_QC_GWAS/check_heterozygosity_rate.R
Rscript --no-save /sc/arion/projects/tauomics/Shrishtee/1_QC_GWAS/heterozygosity_outliers_list.R
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt
plink --bfile PART_KatQC_Hg19.6 --remove het_fail_ind.txt --make-bed --out PART_KatQC_Hg19.7

```

```bash
#Checks for relationships between individuals with a pihat > 0.2.
plink --bfile PART_KatQC_Hg19.7 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2

awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome > zoom_pihat.genome


```

### Prepare data prior to submission to Imputation server


```bash
#LOAD REQUIRED MODULES

ml plink
ml vcftools
ml bcftools #To use bgzip
ml minimac4

dir=/sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/KAT
```

1. Make the .frq file

```bash

plink --bfile ${dir}/PART_KatQC_Hg19.7 --keep-allele-order --freq --out ${dir}/PART_KatQC_Hg19.7 --allow-no-sex
```
2. Separate into chromosome

```bash
for chr in {1..22} X Y; do
    # Extract variants for this chromosome
    plink --bfile ${dir}/PART_KatQC_Hg19.7 --chr $chr --make-bed --out ${dir}/chrom/PART_KatQC_Hg19_separate.$chr
done
```

3. Make vcf files

```bash

for chr in {1..22};
  do
  plink --bfile ${dir}/chrom/PART_KatQC_Hg19_separate.$chr --recode vcf --chr $chr --out ${dir}/chrom/PART_KatQC_Hg19_separate.$chr --allow-no-sex
done
```

4. Sort and bgzip

```bash
for chr in {1..22};
  do
  vcf-sort ${dir}/chrom/PART_KatQC_Hg19_separate.$chr.vcf | bgzip -c > ${dir}/chrom/pre_impute_PART_KatQC_Hg19_separate.$chr.vcf.gz
done
```

5. Index vcf files

```bash
for chr in {1..22};
  do
  tabix -p vcf ${dir}/chrom/pre_impute_PART_KatQC_Hg19_separate.$chr.vcf.gz
done
```

6. Submit to the imputation server or use minimac4/1.6 (latest version)

```bash
msav_files=/sc/arion/projects/tauomics/Shrishtee/Imputation/g1k_p3_msav_files_with_estimates

for chr in {1..22};
  do
  minimac4 ${msav_files}/1000g_phase3_v5.chr$chr.with_parameter_estimates.msav ${dir}/chrom/pre_impute_PART_KatQC_Hg19_separate.$chr.vcf.gz -o ${dir}/Imputed/imputed_PART_KatQC_Hg19_separate.$chr.vcf.gz
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

```bash
ml minimac4

i=${1}

dir=/sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/KAT
msav_files=/sc/arion/projects/tauomics/Shrishtee/Imputation/g1k_p3_msav_files_with_estimates

minimac4 ${msav_files}/1000g_phase3_v5.chr$i.with_parameter_estimates.msav ${dir}/chrom/pre_impute_PART_KatQC_Hg19_separate.$i.vcf.gz -o ${dir}/Imputed/imputed_PART_KatQC_Hg19_separate.$i.vcf.gz

```

# 3. GWAS

Now, we take the Imputed files and run an association analysis using plink.
* Plink version: plink/1.90b6.21

```bash
i=${1}

ml plink

cd WORKING_DIRECTORY

plink \
--bfile PLINK_FILE$i \
--logistic \
--ci 0.95 \
--pheno Covariate_file \
--pheno-name status \
--covar Covariate_file \
--covar-name PC01,PC02,PC03,sex,express,gsa,ucla,mssm,Haplotype \
--keep-allele-order \
--allow-no-sex \
--adjust \
--out Output_file_chr$i

awk 'NR==1 || $5 == "ADD" {print}' Output_file_chr$i > Output_file_chr$i.ADD


cat  Output_file_chr*.ADD | head -n1 >  Output_file_all_chrom.ADD_gwas
cat Output_file_chr*.ADD | grep -v "SE" >> Output_file_all_chrom.ADD_gwas
cat Output_file_all_chrom.ADD_gwas | grep -v "NA" > Output_file_all_chrom.ADD_gwas_without_NA

```

#### Prep the association analysis results to upload to Locus Zoom
```Rscript Association_analysis/Locus_zoom_prep_and_Manhattan_plot.R ```



