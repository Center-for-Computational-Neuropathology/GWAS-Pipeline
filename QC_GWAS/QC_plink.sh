"""
Created on Mon Jan 29 12:22 2024
Goal: Perform SNP and Sample level QC
@@author: Shrishtee Kandoi (github: ShrishteeKandoi)
Username: kandos01
"""

#Load required modules
ml plink

#Read each chrom per batch
i=${1}

########################################
################ Step 1 ################
########################################

#Investigate missingness per individual and per SNP and make histograms
plink --bfile upenn_ucla_mssm.chr$i --missing
# output: plink.imiss and plink.lmiss, these files show respectively the proportion of missing SNPs per individual and the proportion of missing individuals per SNP.

#Generate plots to visualize the missingness results.
Rscript --no-save hist_miss.R

#Delete SNPs with missingness >0.2
plink --bfile upenn_ucla_mssm.chr$i --geno 0.2 --make-bed --out upenn_ucla_mssm.chr$i.2 –keep-allele-order

#Delete individuals with missingness >0.2
plink --bfile upenn_ucla_mssm.chr$i.2 --mind 0.2 --make-bed --out upenn_ucla_mssm.chr$i.3 –keep-allele-order

#Delete SNPs with missingness >0.02
plink --bfile upenn_ucla_mssm.chr$i.3 --geno 0.02 --make-bed --out upenn_ucla_mssm.chr$i.4 –keep-allele-order

#Delete individuals with missingness >0.02
plink --bfile upenn_ucla_mssm.chr$i.4 --mind 0.02 --make-bed --out upenn_ucla_mssm.chr$i.5 –keep-allele-order

########################################
################ Step 2 ################
########################################

# Check for sex discrepancy.
# Subjects who were a priori determined as females must have a F value of <0.2, and subjects who were a priori determined as males must have a F value >0.8. This F value is based on the X chromosome inbreeding (homozygosity) estimate.
# Subjects who do not fulfil these requirements are flagged "PROBLEM" by PLINK.

plink --bfile upenn_ucla_mssm.chr$i.5 --check-sex --allow-no-sex

# Generate plots to visualize the sex-check results.
Rscript --no-save gender_check.R
#Error: --check-sex/--impute-sex requires at least one polymorphic X chromosome locus.
# These checks indicate that there is one woman with a sex discrepancy, F value of 0.99. (When using other datasets often a few discrepancies will be found). 
# The following two scripts can be used to deal with individuals with a sex discrepancy.
# Note, please use one of the two options below to generate the bfile hapmap_r23a_6, this file we will use in the next step of this tutorial.

# 1) Delete individuals with sex discrepancy.
grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt
# This command generates a list of individuals with the status �PROBLEM�.
plink --bfile upenn_ucla_mssm.chr$i.5 --remove sex_discrepancy.txt --make-bed --out upenn_ucla_mssm.chr$i.6 
# This command removes the list of individuals with the status �PROBLEM�.

# 2) Impute-sex.
#plink --bfile HapMap_3_r3_5 --impute-sex --make-bed --out HapMap_3_r3_6
# This imputes the sex based on the genotype information into your data set.

########################################
################ Step 3 ################
########################################

# Generate a bfile with autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF).

# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' upenn_ucla_mssm.chr$i.6 > snp_1_22.txt
plink --bfile upenn_ucla_mssm.chr$i.6 --extract snp_1_22.txt --make-bed --out upenn_ucla_mssm.chr$i.7

# Generate a plot of the MAF distribution.
plink --bfile upenn_ucla_mssm.chr$i.7 --freq --out MAF_check
Rscript --no-save MAF_check.R

# Remove SNPs with a low MAF frequency.
plink --bfile upenn_ucla_mssm.chr$i.7 --maf 0.05 --make-bed --out upenn_ucla_mssm.chr$i.8
# 1073226 SNPs are left
# A conventional MAF threshold for a regular GWAS is between 0.01 or 0.05, depending on sample size.

########################################
################ Step 4 ################
########################################

# Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).
# Check the distribution of HWE p-values of all SNPs.

plink --bfile upenn_ucla_mssm.chr$i.8 --hardy

#Selecting SNPs with HWE p-value below 0.00001, required for one of the two plot generated by the next Rscript, allows to zoom in on strongly deviating SNPs. 
awk '{ if ($9 <0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe
Rscript --no-save hwe.R

#By default the --hwe option in plink only filters for controls.
#Therefore, we use two steps, first we use a stringent HWE threshold for controls, followed by a less stringent threshold for the case data.
plink --bfile upenn_ucla_mssm.chr$i.8 --hwe 1e-6 --make-bed --out upenn_ucla_mssm_hwe.chr$i.8

# The HWE threshold for the cases filters out only SNPs which deviate extremely from HWE. 
# This second HWE step only focuses on cases because in the controls all SNPs with a HWE p-value < hwe 1e-6 were already removed
plink --bfile upenn_ucla_mssm_hwe.chr$i.8 --hwe 1e-10 --hwe-all --make-bed --out upenn_ucla_mssm_hwe.chr$i.9

# Theoretical background for this step is given in our accompanying article: https://www.ncbi.nlm.nih.gov/pubmed/29484742 .

########################################
################ Step 5 ################
########################################

# Generate a plot of the distribution of the heterozygosity rate of your subjects.
# And remove individuals with a heterozygosity rate deviating more than 3 sd from the mean.

# Checks for heterozygosity are performed on a set of SNPs which are not highly correlated.
# Therefore, to generate a list of non-(highly)correlated SNPs, we exclude high inversion regions (inversion.txt [High LD regions]) and prune the SNPs using the command --indep-pairwise�.
# The parameters �50 5 0.2� stand respectively for: the window size, the number of SNPs to shift the window at each step, and the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously.

plink --bfile upenn_ucla_mssm_hwe.chr$i.9 --exclude inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
# Note, don't delete the file indepSNP.prune.in, we will use this file in later steps of the tutorial.

plink --bfile upenn_ucla_mssm_hwe.chr$i.9 --extract indepSNP.prune.in --het --out R_check
# This file contains your pruned data set.

# Plot of the heterozygosity rate distribution
Rscript --no-save check_heterozygosity_rate.R

# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
# For data manipulation we recommend using UNIX. However, when performing statistical calculations R might be more convenient, hence the use of the Rscript for this step:
Rscript --no-save heterozygosity_outliers_list.R

# Output of the command above: fail-het-qc.txt .
# When using our example data/the HapMap data this list contains 2 individuals (i.e., two individuals have a heterozygosity rate deviating more than 3 SD's from the mean).
# Adapt this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns.
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

# Remove heterozygosity rate outliers.
plink --bfile upenn_ucla_mssm_hwe.chr$i.9 --remove het_fail_ind.txt --make-bed --out upenn_ucla_mssm_hwe.chr$i.10

########################################
################ Step 6 ################
########################################

# It is essential to check datasets you analyse for cryptic relatedness.
# Assuming a random population sample we are going to exclude all individuals above the pihat threshold of 0.2 in this tutorial.

# Check for relationships between individuals with a pihat > 0.2.
plink --bfile upenn_ucla_mssm_hwe.chr$i.10 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2

# The HapMap dataset is known to contain parent-offspring relations. 
# The following commands will visualize specifically these parent-offspring relations, using the z values. 
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome > zoom_pihat.genome

# Generate a plot to assess the type of relationship.
Rscript --no-save Relatedness.R

# The generated plots show a considerable amount of related individuals (explentation plot; PO = parent-offspring, UN = unrelated individuals) in the Hapmap data, this is expected since the dataset was constructed as such.
# Normally, family based data should be analyzed using specific family based methods. In this tutorial, for demonstrative purposes, we treat the relatedness as cryptic relatedness in a random population sample.
# In this tutorial, we aim to remove all 'relatedness' from our dataset.
# To demonstrate that the majority of the relatedness was due to parent-offspring we only include founders (individuals without parents in the dataset).
plink --bfile upenn_ucla_mssm_hwe.chr$i.10 --filter-founders --make-bed --out upenn_ucla_mssm_hwe.chr$i.11

# Now we will look again for individuals with a pihat >0.2.
plink --bfile upenn_ucla_mssm_hwe.chr$i.11 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders
# The file 'pihat_min0.2_in_founders.genome' shows that, after exclusion of all non-founders, only 1 individual pair with a pihat greater than 0.2 remains in the HapMap data.
# This is likely to be a full sib or DZ twin pair based on the Z values. Noteworthy, they were not given the same family identity (FID) in the HapMap data.

# For each pair of 'related' individuals with a pihat > 0.2, we recommend to remove the individual with the lowest call rate. 
plink --bfile upenn_ucla_mssm_hwe.chr$i.11 --missing
# Use an UNIX text editor (e.g., vi(m) ) to check which individual has the highest call rate in the 'related pair'. 

# Generate a list of FID and IID of the individual(s) with a Pihat above 0.2, to check who had the lower call rate of the pair.
# In our dataset the individual 13291  NA07045 had the lower call rate.
vi 0.2_low_call_rate_pihat.txt
i 
13291  NA07045
# Press esc on keyboard! :x
# Press enter on keyboard
# In case of multiple 'related' pairs, the list generated above can be extended using the same method as for our lone 'related' pair.

# Delete the individuals with the lowest call rate in 'related' pairs with a pihat > 0.2 
plink --bfile upenn_ucla_mssm_hwe.chr$i.11 --remove 0.2_low_call_rate_pihat.txt --make-bed --out upenn_ucla_mssm_hwe.chr$i.12

################################################################################################################################

# CONGRATULATIONS!! You've just succesfully completed the first tutorial! You are now able to conduct a proper genetic QC. 

# For the next tutorial, using the script: 2_Main_script_MDS.txt, you need the following files:
# - The bfile HapMap_3_r3_12 (i.e., HapMap_3_r3_12.fam,HapMap_3_r3_12.bed, and HapMap_3_r3_12.bim
# - indepSNP.prune.in

