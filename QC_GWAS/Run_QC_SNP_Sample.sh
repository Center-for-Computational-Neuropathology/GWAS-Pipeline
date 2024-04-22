script=/sc/arion/projects/tauomics/Shrishtee/QC_GWAS/QC_SNP_Sample.sh
Exp=QC

for i in {14..14}; do
    echo "
    #!/bin/bash
    #BSUB -P acc_sharpa01a
    #BSUB -n 8
    #BSUB -W 40:00
    #BSUB -R rusage[mem=48000]
    #BSUB -q premium
    #BSUB -J “QC”
    #BSUB -o logs.$Exp/$i.stdout
    #BSUB -eo logs.$Exp/$i.stderr
    mkdir -p logs.$Exp
    rm -f logs.$Exp/$i.stdout
    rm -f logs.$Exp/$i.stdout
    module purge
    sh $script $i" > run.$i
    bsub <  run.$i
    rm run.$i
done

#sh Run_QC_SNP_Sample.sh

############### Run with named arguments ############

#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 --dir <directory> --geno <geno_value> --mind <mind_value> --Exp <experiment_name> --plinkfilename <plink_filename>"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --dir)
            dir="$2"
            shift 2
            ;;
        --geno)
            geno="$2"
            shift 2
            ;;
        --mind)
            mind="$2"
            shift 2
            ;;
        --Exp)
            Exp="$2"
            shift 2
            ;;
        --plinkfilename)
            plinkfilename="$2"
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done

# Check if all required arguments are provided
if [ -z "$dir" ] || [ -z "$geno" ] || [ -z "$mind" ] || [ -z "$Exp" ] || [ -z "$plinkfilename" ]; then
    usage
fi

script="/sc/arion/projects/tauomics/Shrishtee/QC_GWAS/QC_SNP_Sample.sh"

echo "Running QC for dir: $dir, geno: $geno, mind: $mind, Exp: $Exp, plinkfilename: $plinkfilename"

# Execute QC_SNP_Sample.sh script
echo "
#!/bin/bash
#BSUB -P acc_sharpa01a
#BSUB -n 8
#BSUB -W 40:00
#BSUB -R rusage[mem=48000]
#BSUB -q premium
#BSUB -J $Exp
#BSUB -o logs.$Exp/$plinkfilename.stdout
#BSUB -eo logs.$Exp/$plinkfilename.stderr
mkdir -p logs.$Exp
module purge
sh $script --dir $dir --geno $geno --mind $mind --plinkfilename $plinkfilename
" | bsub

# sh Run_QC_SNP_Sample.sh --dir /path/to/dir --geno 0.2 --mind 0.2 --Exp QC
