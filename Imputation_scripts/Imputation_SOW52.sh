#!/bin/bash
#BSUB -P acc_sharpa01a
#BSUB -n 8
#BSUB -W 40:00
#BSUB -R rusage[mem=48000]
#BSUB -R "span[hosts=1]"
#BSUB -q premium
#BSUB -J “imputation”
#BSUB -o imputation_%J.stdout
#BSUB -e imputation_%J.

ml minimac4

i=${1}

dir=/sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/SOW52
msav_files=/sc/arion/projects/tauomics/Shrishtee/Imputation/g1k_p3_msav_files_with_estimates

minimac4 ${msav_files}/1000g_phase3_v5.chr$i.with_parameter_estimates.msav ${dir}/chrom/pre_impute_PART_SOW52_PART.separate.$i.vcf.gz -o ${dir}/Imputed/imputed_PART_SOW52_PART.separate.$i.vcf.gz
