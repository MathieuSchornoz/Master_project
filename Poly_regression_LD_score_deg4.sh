#! /bin/bash

############

# RUN INFO #

###########


#SBATCH --job-name=polynomial_regression 				# Job name
#SBATCH --mail-type=ALL 					# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=schornozmathieu@gmail.com 			# Where to send mail
#SBATCH --ntasks=1 						# Run on a single core
#SBATCH --time=24:00:00 					# Time limit hrs:min:sec, Define how long the job will run
#SBATCH --nodes=1						# Defines number of node required
#SBATCH --account=sgg
#SBATCH --partition=sgg						# Defines the partition on which the job shall run, may be omitted
#SBATCH --mem=30GB						# Memory required per node
#SBATCH --output=Polynomial_regression_degree4_%j.out			# Standard output and error log


##########

# CALCULATE POLYNOMIAL REGRESSION #

##########

# direction where Script and all Tagging files are located
export DIR='/data/sgg/mathieu/LDAK_project/Tagging/'

(cd $DIR ;
Rscript Poly_LD_score_degree4.R
)

echo "Job done"

# Finish the script
exit 0

