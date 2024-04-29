#!/bin/bash
#BSUB -P acc_sharpa01a
#BSUB -n 8
#BSUB -W 10:00
#BSUB -R rusage[mem=30000]
#BSUB -R "span[hosts=1]"
#BSUB -q premium
#BSUB -J “Crossmap”
#BSUB -o Crossmap_%J.stdout
#BSUB -e Crossmap_%J.stderr

#Read each chrom per batch
# conda init
# conda activate mamba

cd /sc/arion/projects/tauomics/Shrishtee/
ml bcftools

i=${1}

path_dir=PART_DATA_SOW_KAT_CHOP/SOW52/Imputed

CrossMap vcf PART_DATA_SOW_KAT_CHOP/hg19ToHg38.over.chain.gz \
    ${path_dir}/imputed_PART_SOW52_PART.separate.${i}.vcf.gz \
    Reference_files/hg38.fa.gz \
    ${path_dir}/hg38/batch/PART_SOW52_Imputed.chr${i}_hg38.vcf

bcftools sort -o ${path_dir}/hg38/batch/PART_SOW52_Imputed.chr${i}_hg38_sorted.vcf ${path_dir}/hg38/batch/PART_SOW52_Imputed.chr${i}_hg38.vcf
bgzip ${path_dir}/hg38/batch/PART_SOW52_Imputed.chr${i}_hg38_sorted.vcf
tabix -p vcf ${path_dir}/hg38/batch/PART_SOW52_Imputed.chr${i}_hg38_sorted.vcf.gz
