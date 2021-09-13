#! /bin/bash

############

# RUN INFO #

###########


#SBATCH --job-name=multivariate_estimation_direct                           # Job name
#SBATCH --mail-type=ALL                                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=schornozmathieu@gmail.com                   # Where to send mail
#SBATCH --ntasks=1                                              # Run on a single core
#SBATCH --time=48:00:00                                         # Time limit hrs:min:sec, Define how long the job will run
#SBATCH --nodes=1                                               # Defines number of node required
#SBATCH --account=sgg
#SBATCH --partition=sgg                                         # Defines the partition on which the job shall run, may be omitted
#SBATCH --mem=30GB                                              # Memory required per node
#SBATCH --output=Log_multivariate_estimation_direct_%j.out                   # Standard output and error log


#################
#   VARIABLES   #
#################

echo "SLURM_JOB_ID: " $SLURM_JOB_ID
echo "SLURM_ARRAY_ID: " $SLURM_ARRAY_ID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID



##########

# MULTIVARIATE ESTIMATION OF SELECTION PARAMETERS #

##########

# Nelder-Mead
# CG
# L-BFGS-B


export DIR='/data/sgg/mathieu/LDAK_project/Multivariate_estimation'
export HOME='/home/ma8758'

# Change direction of output file
cd $DIR

Rscript multivariate_estimation_all_together.R $DIR/Beta/Z_score_table.txt $DIR/Correlation_data/data-2021-06-23.csv $DIR/4080_irnt_gwas_imputed_v3_both_sexes_LD_coefficient.txt CG

echo "Job done"

# Finish the script
exit 0

