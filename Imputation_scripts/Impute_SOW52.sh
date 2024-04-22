#!/bin/bash
#BSUB -P acc_sharpa01a
#BSUB -n 8
#BSUB -W 40:00
#BSUB -R rusage[mem=48000]
#BSUB -R "span[hosts=1]"
#BSUB -q premium
#BSUB -J “imputation”
#BSUB -o imputation_%J.stdout
#BSUB -e imputation_%J.stderr

cd /sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/SOW52/chrom_subset/
dir=/sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/SOW52/Imputed_files

#Load modules
module purge
ml plink
ml vcftools
ml bcftools
ml minimac4/1.6

#Read each chrom per batch
i=${1}

#Making the ".frq" file
plink --bfile PART_SOW52_$i --keep-allele-order --allow-no-sex --freq --out ${dir}/PART_SOW52_$i

#Making the "vcf" file
plink --bfile PART_SOW52_$i --keep-allele-order --allow-no-sex --recode vcf --out ${dir}/PART_SOW52_$i

#Sort the vcf file and "bgzip" it
vcf-sort ${dir}/PART_SOW52_$i.vcf | bgzip -c >  ${dir}/PART_SOW52_$i.vcf.gz
tabix -p vcf ${dir}/PART_SOW52_$i.vcf.gz

#RUN Imputation
msav_files=/sc/arion/projects/tauomics/Shrishtee/Imputation/g1k_p3_msav_files_with_estimates
minimac4 ${msav_files}/1000g_phase3_v5.chr$i.with_parameter_estimates.msav ${dir}/PART_SOW52_$i.vcf.gz -o $dir/Imputed_PART_SOW52_$i.vcf.gz -t 8 --sites $dir/. --temp-prefix $dir/.
plink --vcf ${dir}/Imputed_PART_SOW52_$i.vcf.gz --recode --nonfounders --allow-no-sex --make-bed --out ${dir}/Imputed_PART_SOW52_$i --double-id

