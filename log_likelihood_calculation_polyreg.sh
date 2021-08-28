#! /bin/bash

############

# RUN INFO #

###########


#SBATCH --job-name=Loglikelihood_calculation_deg4 			# Job name
#SBATCH --mail-type=ALL 						# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=schornozmathieu@gmail.com 				# Where to send mail
#SBATCH --ntasks=1 							# Run on a single core
#SBATCH --time=24:00:00 						# Time limit hrs:min:sec, Define how long the job will run
#SBATCH --nodes=1							# Defines number of node required
#SBATCH --account=sgg
#SBATCH --partition=sgg							# Defines the partition on which the job shall run, may be omitted
#SBATCH --mem=10GB							# Memory required per node
#SBATCH --output=Loglikelihood_calculation_polynomial_regre_deg4_%j.out		# Standard output and error log

# !!!!! ADAPT THE NUMBER OF ARRAY TO THE NUMBER OF FILE TO PROCESS

#SBATCH --array=0-53



##########

# VARIABLES

##########

# Absolute path where the script Log_likelihood_calculation_hpc1_polyreg.R is located
export DIR='/data/sgg/mathieu/LDAK_project'
export HOME='/home/ma8758'

# Path where all summary statistics file are located
FILES=($DIR/Summary_Stat/*both_sexes.txt)
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

##########

# CALCULATE LOGLIKELIHOOD WITH LD ESTIMATION #
# of all the summary statistics files

##########

# Change direction of output files
cd $DIR

# Loop around all summary statistic files
echo ${FILE}
Rscript $DIR/Log_likelihood_calculation_hpc1_polyreg.R $DIR/Tagging/unweighted_LD_regression_deg_2/LDAK+Alpha.20_LD_score.tagging ${FILE}

echo "Job done"

# Finish the script
exit 0

