#LOAD REQUIRED MODULES

ml plink
ml vcftools
ml bcftools #To use bgzip
ml minimac4

dir=/sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/KAT

#Make your .frq file
plink --bfile ${dir}/PART_KatQC_Hg19.7 --keep-allele-order --freq --out ${dir}/PART_KatQC_Hg19.7 --allow-no-sex

#Separate into chromosome
for chr in {1..22} X Y; do
    # Extract variants for this chromosome
    plink --bfile ${dir}/PART_KatQC_Hg19.7 --chr $chr --make-bed --out ${dir}/chrom/PART_KatQC_Hg19_separate.$chr
done

#Make vcf files
# plink --bfile ${dir}/PART_KatQC_Hg19.7 --recode vcf --out ${dir}/PART_KatQC_Hg19.7 --keep-allele-order

for chnum in {1,..22};
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



#BY CHROMOSOME

# for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
#   do
#   plink --bfile upenn_ucla_mssm.chr$chnum.5 --keep-allele-order --freq --out upenn_ucla_mssm_impute.chr$chnum --allow-no-sex
# done

# for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
#   do
#   plink --bfile YOURFILE-updated-chr$chnum --recode vcf --chr $chnum --out YOURFILE$chnum --allow-no-sex
# done

# for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
#   do
# 	vcf-sort YOURFILE$chnum.vcf | bgzip -c >  pre_impute_YOURFILE_$chnum.vcf.gz
# done


#SOW52

dir=/sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/SOW52

#Make your .frq file
plink --bfile ${dir}/SOW52_PART_sex.update.7 --keep-allele-order --freq --out ${dir}/SOW52_PART_sex.update.7 --allow-no-sex

#Separate into chromosome
for chr in {1..22} X Y; do
    # Extract variants for this chromosome
    plink --bfile ${dir}/SOW52_PART_sex.update.7 --chr $chr --make-bed --out ${dir}/chrom/SOW52_PART.separate.$chr
done

#Make vcf files
# plink --bfile ${dir}/PART_KatQC_Hg19.7 --recode vcf --out ${dir}/PART_KatQC_Hg19.7 --keep-allele-order

for chr in {1..22};
  do
  plink --bfile ${dir}/chrom/SOW52_PART.separate.$chr --recode vcf --chr $chr --out ${dir}/chrom/SOW52_PART.separate.$chr --allow-no-sex
done

#Then sort and zip
# vcf-sort ${dir}/PART_KatQC_Hg19.7.vcf | bgzip -c > pre_impute_PART_KatQC_Hg19.7.vcf.gz

#BGZIP vcf files
for chr in {1..22};
  do
  vcf-sort ${dir}/chrom/SOW52_PART.separate.$chr.vcf | bgzip -c > pre_impute_PART_SOW52_PART.separate.$chr.vcf.gz
done

#Index vcf files
for chr in {1..22};
  do
  tabix -p vcf ${dir}/chrom/pre_impute_PART_SOW52_PART.separate.$chr.vcf.gz
done