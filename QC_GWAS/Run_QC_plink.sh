script=/sc/arion/projects/tauomics/Shrishtee/QC_GWAS/QC_plink.sh
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

#sh run_imputation_batch_jobs.sh