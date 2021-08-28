
############################################################

# Log likelihood calculation with 3 parameters and a 4th degree polynomial regression

############################################################

# Install package and library
#install.packages("tydiverse", repos = 'https://stat.ethz.ch/CRAN/')
#install.packages("dplyr", repos = 'https://stat.ethz.ch/CRAN/')
#install.packages("stringr", repos = 'https://stat.ethz.ch/CRAN/')

# Installing packages
x = c("dplyr", "tibble", "stringr")
lapply(x, require, character.only = TRUE)

# Retrieve script input
args = commandArgs(TRUE)
input1 = args[1]
input2 = args[2]

# Open input file as table
Tagging_file = read.table(file = args[1], header = TRUE, sep = " ", stringsAsFactors = FALSE)
Sum_stat = read.table(file = args[2], header = TRUE, sep = " ", stringsAsFactors = FALSE)

# Filtering out last row of Tagging file
Tagging_file = Tagging_file[1:(nrow(Tagging_file)-2),]

cat("Binding Summary statistics data with Tagging file:\n", input1, "\n")
cat("Summary Statistic:\n", input2, "\n")

####### Binding Summary statistics file with the tagging file
geno_data = Sum_stat %>%
  inner_join(Tagging_file, by = c("Predictor")) %>%
  mutate_at(c("Z", "Tagging", "Weight", "MAF", "Exp_Heritability", "Base"), as.numeric) %>%
  mutate_at(c("Neighbours", "Categories"), as.integer) %>%
  mutate_at(c("Predictor"), as.character)

rm(Sum_stat, Tagging_file)

# Inconsistent condition and filtration
cond1 = (geno_data$A1.x == "A" | geno_data$A2.x == "A") & !(geno_data$A1.y == "A" | geno_data$A2.y == "A")
cond2 = (geno_data$A1.x == "C" | geno_data$A2.x == "C") & !(geno_data$A1.y == "C" | geno_data$A2.y == "C")
cond3 = (geno_data$A1.x == "G" | geno_data$A2.x == "G") & !(geno_data$A1.y == "G" | geno_data$A2.y == "G")
cond4 = (geno_data$A1.x == "T" | geno_data$A2.x == "T") & !(geno_data$A1.y == "T" | geno_data$A2.y == "T")
condition_inconsistent = cond1|cond2|cond3|cond4
rm(cond1,cond2,cond3,cond4)

cat("\n Number of predictors with inconsistent alleles:",nrow(geno_data[condition_inconsistent,]),"\n")
geno_data = geno_data[!condition_inconsistent,] 
rm(condition_inconsistent)

# Ambiguous Condition and filtration
cond1 = (geno_data$A1.x == "A" & geno_data$A2.x == "T") | (geno_data$A1.y == "A" & geno_data$A2.y == "T")
cond2 = (geno_data$A1.x == "C" & geno_data$A2.x == "G") | (geno_data$A1.y == "C" & geno_data$A2.y == "G")
cond3 = (geno_data$A1.x == "G" & geno_data$A2.x == "C") | (geno_data$A1.y == "G" & geno_data$A2.y == "C")
cond4 = (geno_data$A1.x == "T" & geno_data$A2.x == "A") | (geno_data$A1.y == "T" & geno_data$A2.y == "A")
condition_ambiguous = cond1|cond2|cond3|cond4
rm(cond1,cond2,cond3,cond4)

cat("\n Number of predictore with ambiguous alleles:", nrow(geno_data[condition_ambiguous,]),"\n")
#geno_data = geno_data[!condition_ambiguous,]
rm(condition_ambiguous)

# Change tagging value to 1 if less than 1
cat("\n Number of row with 'Tagging' value below 1: ", nrow(geno_data[geno_data$Tagging < 1,]),"\n")
geno_data$Tagging[geno_data$Tagging < 1] = 1

# Further data filtering
# cat(" \n Number of row with negative 'Base' value: ", nrow(geno_data[geno_data$Base <= 0,]),"\n")
#  geno_data = geno_data %>%
#    filter(Base > 0, Z != 0)

cat("\n Structure of geno_data: \n")
str(geno_data)
summary(geno_data)

# Check degree of polynomial regression
# and calculation of Base according to S value
if (names(geno_data[ncol(geno_data)]) == "LD_quadratic" ) {
  print("Estimation of LD score follows a quadratic polynomial regression")
  
  Base_calculation = function(alpha_parameter){
    Base_vector = geno_data$LD_Intercept + geno_data$LD_linear * alpha_parameter + geno_data$LD_quadratic * alpha_parameter^2
    return(Base_vector)
  }
  
} else if (names(geno_data[ncol(geno_data)]) == "LD_quartic") {
  print("Estimation of LD score follows a X^4 polynomial regression")
  
  Base_calculation = function(alpha_parameter){
    Base_vector = geno_data$LD_Intercept + geno_data$LD_linear * alpha_parameter + geno_data$LD_quadratic * alpha_parameter^2 +geno_data$LD_cubic * alpha_parameter^3 + geno_data$LD_quartic * alpha_parameter^4
    return(Base_vector)
  }
} else {
  print("Cannot recognize degree of polynomial regression lf LD score")
}

#Test of Base calculation
test = Base_calculation(0.1)
head(test)
head(geno_data$Base)

# Definition of Log_likelihood function
log_likelihood = function(par, table_of_interest){
  # Parameters
  A = par[1]
  tau = abs(par[2])
  alpha = par[3]
  my_data_table = data.frame(table_of_interest)
  # Sample size
  n = my_data_table[1,"n"]
  # Calculate E[J_s]
  #Expect_S = A + n * my_data_table[colnames(my_data_table) == "Base"] * abs(tau)
  Expect_S = A + n * Base_calculation(alpha) * abs(tau)
  # Inverted LD score
  inv_LD_score = 1/my_data_table[colnames(my_data_table) == "Tagging"]
  # Z-squared
  S = (my_data_table[colnames(my_data_table) == "Z"])^2
  # Compute log-likelihood function
  loglss = inv_LD_score * (-S/(2*Expect_S) -log(S)/2 - log(2*Expect_S)/2 - log(2*pi)/2)
  
  cat(sum(is.na(loglss)), ",", nrow(loglss), ",", A, ",", tau, ",",alpha, ",", sum(loglss, na.rm = F), "\n")
  return(sapply(loglss, mean, na.rm = T) * nrow(loglss)) # or length
}

################################################
# Optimization of Log_likelihood function 

# Define parameter for optimization
ratio = c(1,1e-5,1)

# Optimization using A=1 and tau=0 as starting point
start_time = Sys.time()
optim_result = optim(par = c(1,0,0), fn = log_likelihood, table_of_interest = geno_data,
                     control = list(parscale = ratio, fnscale = -1))
end_time = Sys.time()
end_time - start_time
# Control of optimization
#result_log_call = log_likelihood(c(optim_result$par), geno_data)

# Print optimized parameters value
cat("\n ******************************************** \n", 
    format("Optimized log_likelihood function results: ", justify = "left", width = 14),
    "\n  Parameters A, Tau and S: ", format(optim_result$par, width = 8), 
    "\n  Value of log_likelihood: ", format(optim_result$value, width = 8), 
    "\n  Convergence (0 is yes): ", format(optim_result$convergence, width = 8),
    "\n ******************************************** \n")

# Make datafram from the list output of optim
x = gsub("[A-z \\.\\(\\)]", "", input1)
list_names = c("Intercept", "Tau", "Alpha", "value", "counts_function", "count_gradient", "convergence")
optim_result_frame = data.frame(matrix(unlist(optim_result), nrow = 1, byrow = TRUE), stringsAsFactors=FALSE)
colnames(optim_result_frame) = list_names
optim_result_frame$Tau = abs(optim_result_frame$Tau)
optim_result_frame = optim_result_frame %>%
  mutate( Sample_size = geno_data[1,"n"],
          Heritability = abs(optim_result_frame$Tau) * sum(geno_data$Exp_Heritability))

# Prepare filename for exportation
input1 = unlist(strsplit(input1, "/"))
input1 = input1[length(input1)]
input2 = unlist(strsplit(input2, "/"))
input2 = input2[length(input2)]
input2 = str_replace(input2, ".txt", "")
input2 = gsub('\\.', '_', input2)
x = c("Log_lik", input2, input1, ".txt")
filename = paste(x, collapse = "_")
filename = str_replace(filename, ".tagging_", "")
filename = gsub('imputed_v3_both_sexes_LDAK_','', filename)

# Prepare a column containing the phenotype name according genetic correlation file from Neale Lab
phenotype_name = input2
phenotype_name = gsub('_gwas_imputed_v3_both_sexes','',phenotype_name)

# Check for the presence of "irnt" sequence in the name and erase it
phenotype_name = str_remove(phenotype_name, "_irnt")

# Integrate the phenotype name in the output file
optim_result_frame = optim_result_frame %>%
  mutate(Phenotype = phenotype_name)

cat("\n Filename is:", filename, "\n")

# Do I have to create a table where the regression coefficient for all SNP are saved ?
# Let's do this it might be useful in the future
# Create a table with the coefficient of the LD regression

LD_regression_coefficient = geno_data %>%
  select(Predictor,Z, LD_Intercept, LD_linear, LD_quadratic, LD_cubic, LD_quartic) 
filename_LD = c(input2, "LD_coefficient.txt")
filename_LD = paste(filename_LD, collapse = "_")

# Creation of table
write.table(optim_result_frame, file = filename, sep = ",",  row.names = FALSE, )
write.table(LD_regression_coefficient, file = filename_LD, sep = ",", row.names = FALSE)

cat("Done for file: ", input2, "\n")
print("****************************************************************************************")
