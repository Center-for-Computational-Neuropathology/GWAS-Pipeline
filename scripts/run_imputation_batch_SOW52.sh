script=/sc/arion/projects/tauomics/Shrishtee/Imputation/scripts/Imputation_SOW52.sh
Exp=imputation

for i in {1..22}; do
    echo "
    #!/bin/bash
    #BSUB -P acc_tauomics
    #BSUB -n 8
    #BSUB -W 40:00
    #BSUB -R rusage[mem=48000]
    #BSUB -q premium
    #BSUB -J “Imputation”
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

#sh run_imputation_batch_jobs.sh