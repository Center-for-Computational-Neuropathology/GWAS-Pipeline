#!/bin/bash
#BSUB -P PROJECT_ID
#BSUB -n 8
#BSUB -W 40:00
#BSUB -R rusage[mem=32000]
#BSUB -q premium
#BSUB -J "Somalier"
#BSUB -o Somalier%J.stdout
#BSUB -e Somalier%J.stderr

# =========================================================
# Somalier Ancestry and Relatedness Pipeline
# =========================================================

# Load required modules
ml bcftools
ml plink2
ml htslib

# Set working directory
WORKDIR="$PWD"
mkdir -p ${WORKDIR}/VCF_Files
mkdir -p ${WORKDIR}/somalier_output/extracted

# Set file paths and references
FILENAME="your_genotype_file"  # Base filename for your genotype data
SITES="path/to/somalier/sites.hg19.vcf.gz"  # Somalier sites file
REF="path/to/reference/hg19.fa.gz"  # Reference genome
ANCESTRY_LABELS="path/to/somalier/ancestry-labels-1kg.tsv"  # Ancestry labels file
LABELED_SAMPLES="path/to/somalier/1kg-somalier"  # 1000 Genomes somalier files
VCF_DIR="${WORKDIR}/VCF_Files"
SOMALIER_OUT="${WORKDIR}/somalier_output"
EXTRACTED="${SOMALIER_OUT}/extracted"

# =========================================================
# STEP 1: Convert PLINK to VCF
# =========================================================
echo "Converting PLINK files to VCF format..."

# Create VCF file from PLINK binary files
plink2 --bfile ${WORKDIR}/${FILENAME} \
       --ref-from-fa force \
       --fa ${REF} \
       --recode vcf \
       --not-chr X,Y \
       --out ${VCF_DIR}/${FILENAME}

cd ${VCF_DIR}

# Compress and index VCF
bgzip ${FILENAME}.vcf
tabix -p vcf ${FILENAME}.vcf.gz

# Check if chromosomes are named as "1" or "chr1"
CHR_FORMAT=$(bcftools view ${VCF_DIR}/${FILENAME}.vcf.gz | grep -v '^#' | head -n 1 | cut -f1)
echo "Detected chromosome format: ${CHR_FORMAT}"

# If needed, rename chromosomes to match reference format
if [[ ${CHR_FORMAT} =~ ^[0-9]+$ ]]; then
    echo "Renaming chromosomes to chr format..."
    
    # Create chromosome name conversion file
    for i in {1..22}; do
        echo -e "${i}\tchr${i}" >> chr_name_conv.txt
    done
    
    # Rename chromosomes
    bcftools annotate --rename-chrs chr_name_conv.txt ${FILENAME}.vcf.gz -Oz -o ${FILENAME}_renamed.vcf.gz
    tabix -p vcf ${FILENAME}_renamed.vcf.gz
    
    # Set the correctly formatted VCF for somalier
    FINAL_VCF="${FILENAME}_renamed.vcf.gz"
else
    FINAL_VCF="${FILENAME}.vcf.gz"
fi

# =========================================================
# STEP 2: Extract Sites for Somalier
# =========================================================
echo "Extracting sites for Somalier analysis..."

cd ${SOMALIER_OUT}

# Check if somalier executable exists, if not download it
if [ ! -f ./somalier ]; then
    echo "Downloading somalier..."
    wget https://github.com/brentp/somalier/releases/download/v0.2.19/somalier
    chmod +x somalier
fi

# Extract sites needed for somalier
./somalier extract -d ${EXTRACTED} --sites ${SITES} -f ${REF} ${VCF_DIR}/${FINAL_VCF}

# =========================================================
# STEP 3: Run Relatedness Analysis
# =========================================================
echo "Running relatedness analysis..."

cd ${SOMALIER_OUT}
./somalier relate ${EXTRACTED}/*.somalier

# =========================================================
# STEP 4: Run Ancestry Analysis
# =========================================================
echo "Running ancestry analysis..."

# Check file count to determine if batching is needed
FILE_COUNT=$(ls ${EXTRACTED}/*.somalier | wc -l)
echo "Total somalier files: ${FILE_COUNT}"

# If file count is large, split into batches
if [ ${FILE_COUNT} -gt 5000 ]; then
    echo "Large number of files detected. Running ancestry analysis in batches..."
    
    # Create batch directories
    BATCH_SIZE=5000
    BATCH_COUNT=$(( (FILE_COUNT + BATCH_SIZE - 1) / BATCH_SIZE ))
    
    for i in $(seq 1 ${BATCH_COUNT}); do
        mkdir -p "${EXTRACTED}/batch_${i}"
    done
    
    # Get list of all somalier files and distribute to batches
    cd ${EXTRACTED}
    FILES=(*.somalier)
    CURRENT_BATCH=1
    
    for ((i=0; i<${FILE_COUNT}; i++)); do
        if [ $((i % BATCH_SIZE)) -eq 0 ] && [ $i -ne 0 ]; then
            ((CURRENT_BATCH++))
        fi
        cp "${FILES[$i]}" "batch_${CURRENT_BATCH}/"
    done
    
    # Process each batch
    cd ${SOMALIER_OUT}
    for i in $(seq 1 ${BATCH_COUNT}); do
        echo "Processing batch ${i}..."
        BATCH_DIR="${EXTRACTED}/batch_${i}"
        
        ./somalier ancestry \
            --labels "${ANCESTRY_LABELS}" \
            "${LABELED_SAMPLES}"/*somalier \
            ++ \
            "${BATCH_DIR}"/*somalier
        
        # Rename output files with batch number
        if [ -f "somalier-ancestry.somalier-ancestry.html" ] && [ -f "somalier-ancestry.somalier-ancestry.tsv" ]; then
            mv somalier-ancestry.somalier-ancestry.html "somalier-ancestry_batch_${i}.html"
            mv somalier-ancestry.somalier-ancestry.tsv "somalier-ancestry_batch_${i}.tsv"
            echo "Successfully saved results for batch ${i}"
        else
            echo "Warning: Output files not found for batch ${i}"
        fi
    done
else
    # Run ancestry analysis on all files at once
    ./somalier ancestry \
        --labels "${ANCESTRY_LABELS}" \
        "${LABELED_SAMPLES}"/*somalier \
        ++ \
        "${EXTRACTED}"/*.somalier
    
    echo "Ancestry analysis completed."
fi

echo "Pipeline completed successfully."
