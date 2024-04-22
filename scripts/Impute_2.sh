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

cd /sc/arion/projects/tauomics/Shrishtee/PART_CHOP_DATA_23Feb2024/
dir=/sc/arion/projects/tauomics/Shrishtee/PART_CHOP_DATA_23Feb2024/Imputed_files

#Load modules
module purge
ml plink
ml vcftools
ml bcftools
ml minimac4/1.6

#Read each chrom per batch
i=${1}

#Making the ".frq" file
plink --bfile IKF_Neurology_x600_2-22-24_PLINK_FILES.txt.GSA-24v3.ped.10 --keep-allele-order --allow-no-sex --freq --out IKF_Neurology_x600_2-22-24_PLINK_FILES.txt.GSA-24v3.ped.10_frq

#Making the "vcf" file
plink --bfile IKF_Neurology_x600_2-22-24_PLINK_FILES.txt.GSA-24v3.ped.10 --keep-allele-order --allow-no-sex --recode vcf --out IKF_Neurology_x600_2-22-24_PLINK_FILES.txt.GSA-24v3.ped.10_vcf

#Sort the vcf file and "bgzip" it
vcf-sort IKF_Neurology_x600_2-22-24_PLINK_FILES.txt.GSA-24v3.ped.10_vcf.vcf | bgzip -c >  IKF_Neurology_x600_2-22-24_PLINK_FILES.txt.GSA-24v3.ped.10_vcf.vcf.gz
tabix -p vcf IKF_Neurology_x600_2-22-24_PLINK_FILES.txt.GSA-24v3.ped.10_vcf.vcf.gz

#RUN Imputation
minimac4 1000g_phase3_v5.chr$i.with_parameter_estimates.msav IKF_Neurology_x600_2-22-24_PLINK_FILES.txt.GSA-24v3.ped.10_vcf.vcf.gz -o $dir/IKF_Neurology_x600_2-22-24_PLINK_FILES.txt.GSA-24v3.ped.10_Imputed.vcf.gz -t 8 --sites $dir/. --temp-prefix $dir/.