#! /bin/bash

############

# RUN INFO #

###########


#SBATCH --job-name=tagging_calculation 				# Job name
#SBATCH --mail-type=ALL 					# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=schornozmathieu@gmail.com 			# Where to send mail
#SBATCH --ntasks=1 						# Run on a single core
#SBATCH --time=24:00:00 					# Time limit hrs:min:sec, Define how long the job will run
#SBATCH --nodes=1						# Defines number of node required
#SBATCH --account=sgg
#SBATCH --partition=sgg						# Defines the partition on which the job shall run, may be omitted
#SBATCH --mem=20GB						# Memory required per node
#SBATCH --output=/home/ma8758/Projects/LDAK_project/tagging_calculation_BMI_%j.out			# Standard output and error log
#SBATCH --array=0-30


#################
#   VARIABLES   #
#################
​
echo "SLURM_JOB_ID: " $SLURM_JOB_ID
echo "SLURM_ARRAY_ID: " $SLURM_ARRAY_ID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID



##########

# CALCULATE TAGGINGS #

##########

# Produce unweighted grid of tagging file with fixed alpha value (from -1 to 0.5)


alpha=`echo $SLURM_ARRAY_TASK_ID | awk '{print ($1-20)/20}'`; echo $alpha
./ldak5.1.linux --calc-tagging Tagging/LDAK+Alpha_$SLURM_ARRAY_TASK_ID --bfile Reference/ref.BMI --ignore-weights YES --power $alpha --window-kb 1000


echo "Job done for LDAK_BMI_Alpha_" $alpha

# Finish the script
exit 0

