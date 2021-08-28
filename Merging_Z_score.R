###############################################################

# COMBINE Z-STATISTICS FROM ALL LD SCORE FILE

###############################################################


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


cat("Binding Summary statistics:\n", input1, "\n")

# Prepare a column containing the phenotype name according genetic correlation file from Neale Lab
phenotype_name_1 = input1
phenotype_name_1 = gsub('_gwas_imputed_v3_both_sexes_LD_coefficient.txt','',phenotype_name_1)
phenotype_name_1 = gsub('_irnt', '', phenotype_name_1)

First_table = read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)

First_table = First_table[, c("Predictor","Z")]
colnames(First_table) = c("Predictor", phenotype_name_1)
as.numeric(First_table[,2])

Z_score_dataframe = First_table

for (i in 2:length(args) ) {
  
  input2 = args[i]
  phenotype_name_i = input2
  phenotype_name_i = gsub('_gwas_imputed_v3_both_sexes_LD_coefficient.txt','', phenotype_name_i)
  phenotype_name_i = gsub('_irnt', '', phenotype_name_i)
  
  cat("\nStarting with file:", phenotype_name_i, "for the ", i, "iteration \n")
  
  Second_table = read.table(file = args[i], header = TRUE, sep = ",", stringsAsFactors = FALSE)
  Second_table = Second_table[,c("Predictor","Z")]
  colnames(Second_table) = c("Predictor", phenotype_name_i)
  
  #binding first table with input i
  Z_score_dataframe = Z_score_dataframe %>%
    inner_join(Second_table, by = c("Predictor")) %>%
    mutate_at(c(phenotype_name_i), as.numeric)
  
  print(head(Z_score_dataframe))
  
}

cat("\nDone with the iteration \n")

# Export the Z_score file
write.table(Z_score_dataframe, file = "Z_score_table.txt", sep = ",", row.names = FALSE)





