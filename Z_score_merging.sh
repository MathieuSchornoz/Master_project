#! /bin/bash

############

# RUN INFO #

###########


#SBATCH --job-name=Z_score_merging 				# Job name
#SBATCH --mail-type=ALL 					# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=schornozmathieu@gmail.com 			# Where to send mail
#SBATCH --ntasks=1 						# Run on a single core
#SBATCH --time=24:00:00 					# Time limit hrs:min:sec, Define how long the job will run
#SBATCH --nodes=1						# Defines number of node required
#SBATCH --account=sgg
#SBATCH --partition=sgg						# Defines the partition on which the job shall run, may be omitted
#SBATCH --mem=40GB						# Memory required per node
#SBATCH --output=Z_score_output.out				# Standard output and error log


##########

# Z SCORE MERGING #

##########

FILE=$(ls *LD_coefficient.txt)
echo $FILE

Rscript Merging_Z_score.R $FILE


echo "Job done"

# Finish the script
exit 0

