
######################################################

# CREATION FOUR ORDER POLYNOMIAL REGRESSION LD SCORE 

######################################################

# Install package and library
#install.packages("tydiverse", repos = 'https://stat.ethz.ch/CRAN/')
#install.packages("dplyr", repos = 'https://stat.ethz.ch/CRAN/')
#install.packages("stringr", repos = 'https://stat.ethz.ch/CRAN/')

# Installing packages
x = c("dplyr", "tibble", "stringr")
lapply(x, require, character.only = TRUE)

# Retrieve script input
temp = list.files(pattern = "LDAK_weigh")
Tag_list = lapply(temp, read.table, header = TRUE, sep = " ", stringsAsFactors = FALSE)

# Rename each element of List
names(Tag_list) = temp
names(Tag_list)

# Rename each Base column according to S value
for (i in 1:length(Tag_list)) {
  names(Tag_list[i])
  print(names(Tag_list[i]))
  S_value = gsub("LDAK_weighted_Alpha.", "", names(Tag_list[i]))
  S_value = gsub(".tagging", "", S_value)
  print(S_value)
  print(names(Tag_list[i]$Base))
  new_name = c("Base", S_value)
  new_name = paste(new_name, collapse = "_")
  print(new_name)
  
  names(Tag_list[[i]])
  colnames(Tag_list[[i]]) = c("Predictor", "A1", "A2", "Neighbours", "Tagging", 
                              "Weight", "MAF", "Categories", "Exp_Heritability", new_name)
  #print(colnames(Tag_list[[i]]))
}

# Access Base of each SNP via all S value and change the column name according to the S value
LD_matrix = matrix(ncol = length(Tag_list), nrow = length(Tag_list[[1]]$Predictor))
name_col = vector(mode = "expression", length = length(Tag_list))

for (i in 1:length(Tag_list)) {
  LD_matrix[,i] = Tag_list[[i]]$Base
  name_col[i] = colnames(Tag_list[[i]][ncol(Tag_list[[i]])])
  #print(name_col)
  colnames(LD_matrix) = name_col
}

# reorder column of the matrix
new_order = str_sort(colnames(LD_matrix), numeric = TRUE)
LD_matrix = LD_matrix[,new_order]

# small check
sample_LD_matrix = c(504869,4278173,7013357,268264,2411145)
test_LD = LD_matrix[sample_LD_matrix, ]
test_LD
jpeg("LD_score_sample_degree_4.jpg")
matplot(t(test_LD), type = "b", pch = 1, col = 1:nrow(test_LD))
legend("topright", legend = 1:nrow(test_LD), col = 1:nrow(test_LD), pch = 1)
dev.off()

# Apply polynomial regression to each row
point_S = 0:(ncol(LD_matrix)-1)/20 -1
poly_regression = lm(LD_matrix[1,] ~ poly(point_S, degree = 4, raw = TRUE))
summary(poly_regression)

poly_coef = poly_regression$coefficients
poly_coef
names(poly_coef) = c("LD_Intercept", "LD_linear", "LD_quadratic", "LD_cubic", "LD_quartic")

# Apply poly regression to each row
matrix_poly_coef = matrix(ncol = 5, nrow = nrow(LD_matrix))

for (i in 1:nrow(LD_matrix)) {
  poly_regression = lm(LD_matrix[i,] ~ poly(point_S, degree = 4, raw = TRUE))
  poly_coef = poly_regression$coefficients
  names(poly_coef) = c("LD_Intercept", "LD_linear", "LD_quadratic", "LD_cubic", "LD_quartic")
  matrix_poly_coef[i,] = poly_coef
}

# Check size of matrix of coefficient
cat("\n Number of rows of the matrix:", nrow(matrix_poly_coef), "\n")

# Return the object
new_name_col = c("Base", "LD_Intercept", "LD_linear", "LD_quadratic", "LD_cubic", "LD_quartic")
for (i in 1:length(Tag_list)) {
  Tag_list[[i]] = cbind(Tag_list[[i]], matrix_poly_coef, deparse.level = 0)
  names(Tag_list[[i]])[10:15] = new_name_col
  filename = str_replace(names(Tag_list[i]), ".tagging", "_LD_score_degree_4.tagging")
  write.table(Tag_list[[i]], file = filename, sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

