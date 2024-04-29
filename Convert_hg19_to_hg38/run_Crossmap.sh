script=/sc/arion/projects/tauomics/Shrishtee/GWAS_Pipeline_CraryLab/Convert_hg19_to_hg38/Crossmap.sh
Exp=Crossmap

for i in {22..22}; do
    echo "
    #!/bin/bash
    #BSUB -P acc_tauomics
    #BSUB -n 8
    #BSUB -W 10:00
    #BSUB -R rusage[mem=30000]
    #BSUB -q premium
    #BSUB -J "Crossmap"
    #BSUB -o logs.$Exp/chr$i.stdout
    #BSUB -eo logs.$Exp/chr$i.stderr
    mkdir -p logs.$Exp
    rm -f logs.$Exp/$i.stdout
    rm -f logs.$Exp/$i.stdout
    module purge
    sh $script $i" > run.$i
    bsub <  run.$i
    rm run.$i
done
