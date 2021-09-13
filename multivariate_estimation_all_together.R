####################################

# CREATE DATAFRAME CONTAINING TAU, ALPHA, SAMPLE SIZE OF ALL PHENOTYPE

####################################

# Install package and library
# install.packages("tydiverse", repos = 'https://stat.ethz.ch/CRAN/')
# install.packages("dplyr", repos = 'https://stat.ethz.ch/CRAN/')
# install.packages("stringr", repos = 'https://stat.ethz.ch/CRAN/')
# install.packages("corrplot",repos = 'https://stat.ethz.ch/CRAN/')
# install.packages("factoextra",repos = 'https://stat.ethz.ch/CRAN/')
# install.packages(CulsterR,repos = 'https://stat.ethz.ch/CRAN/')
#install.packages("MASS", repos = 'https://stat.ethz.ch/CRAN/')

# load library
#library("dplyr")
#library("readr")

x = c("dplyr", "tibble", "stringr","readr")
lapply(x, require, character.only = TRUE)


# Retrieve script input
args <- commandArgs(trailingOnly = TRUE)

input1 = args[1]
#Z_score
input2 = args[2]
#data correlation
input3 = args[3]  
#LD coefficient
input4 = args[4]

##### Set location of file in server
path_file = "../Results_weighted_deg_4/"
filelist = list.files(path = path_file, pattern = "Log_lik_", )
setwd(path_file)

#assuming CSV with a header    
datalist = lapply(filelist, FUN=read.table, header=TRUE, sep = ",")

#assuming the same header/columns for all files
output_Loglik = do.call("rbind", datalist) 

# Reorder column 
output_Loglik = output_Loglik[,c(10,1:4,8,9,5:7)]

# erase unused file
#rm(list=setdiff(ls(), "output_Loglik"))

####################################

# RETRIEVING GENETIC CORRELATION 

####################################

library(corrplot)
#library(factoextra)
#library(RColorBrewer)

#theme_set(theme_bw())

# get correlation data from linux 
genetic_correlation = read_csv(file = args[2], 
                               col_types = cols(ID1 = col_character()))


## If P-value above 0.05 change rg to zero
above_05 = genetic_correlation$p > 0.05
genetic_correlation[above_05,"rg"] = 0
nrow(genetic_correlation[genetic_correlation$rg == 0,])

## create the matrix of correlation
matrix_correlation = matrix(NA, nrow = 54, ncol = 54)

# matrix with dimension 54x54
colnames(matrix_correlation) = unique(genetic_correlation$ID1)
rownames(matrix_correlation) = unique(genetic_correlation$ID1)

# Remplace all diagonal element by one
diag(matrix_correlation) = 1

# Return correlation matrix
for (iteration in 1:ncol(matrix_correlation)) {
  element_to_paste = seq(from = (iteration-1) * 53 + 1, to = (iteration-1) * 53 + 53, by = 1)
  discard_element = iteration - 1
  
  value_to_paste = as.vector(unlist(genetic_correlation[element_to_paste, "rg"]))
  
  x = c(rep(999, discard_element),1, tail(value_to_paste, nrow(matrix_correlation)-iteration))
  
  matrix_correlation[,iteration] = x
}

# Get symmetrical matrix
matrix_correlation[upper.tri(matrix_correlation)] = t(matrix_correlation)[upper.tri(matrix_correlation)]

head(matrix_correlation)

# Check if matrix of correlation is symmetrical
cat("\nThe matrix is symmetrical ? -", isSymmetric(matrix_correlation), "\n") 

###########################

# Return significance matrix
significance_matrix = matrix(NA, nrow = 54, ncol = 54)
colnames(significance_matrix) = unique(genetic_correlation$ID1)
rownames(significance_matrix) = unique(genetic_correlation$ID1)

#sequence(nvec = ncol(significance_matrix)-1, from = (iteration-1) * 53 + 1, by = 1)

# return significance matrix
for (iteration in 1:ncol(significance_matrix)) {
  element_to_paste = seq(from = (iteration-1) * 53 + 1, to = (iteration-1) * 53 + 53, by = 1)
  discard_element = iteration - 1
  
  value_to_paste = as.vector(unlist(genetic_correlation[element_to_paste, "p"]))
  
  x = c(rep(999, discard_element),0, tail(value_to_paste, nrow(significance_matrix)-iteration))
  
  significance_matrix[,iteration] = x
}

significance_matrix[upper.tri(significance_matrix)] = t(significance_matrix)[upper.tri(significance_matrix)]


####################################

# CREATION OF LAMBDA MATRIX 

####################################

### First, create the Lambda zero matrix

## create the matrix of correlation
size_matrix = ceiling(sqrt(nrow(genetic_correlation)))
Lambda_matrix = matrix(NA, nrow = size_matrix, ncol = size_matrix)

# matrix with dimension 54x54
colnames(Lambda_matrix) = unique(genetic_correlation$ID1)
rownames(Lambda_matrix) = unique(genetic_correlation$ID1)


# Remplace all diagonal element by one
diag(Lambda_matrix) = 1


#sequence(nvec = ncol(Lambda_matrix)-1, from = (iteration-1) * 53 + 1, by = 1)
# Return Lambda matrix
for (iteration in 1:ncol(Lambda_matrix)) {
  element_to_paste = seq(from = (iteration-1) * 53 + 1, to = (iteration-1) * 53 + 53, by = 1)
  discard_element = iteration - 1
  
  value_to_paste = as.vector(unlist(genetic_correlation[element_to_paste, "rg intercept"]))
  
  x = c(rep(999, discard_element),1, tail(value_to_paste, nrow(Lambda_matrix)-iteration))
  
  Lambda_matrix[,iteration] = x
}

# Get symmetrical Lambda_matrix
Lambda_matrix[upper.tri(Lambda_matrix)] = t(Lambda_matrix)[upper.tri(Lambda_matrix)]
head(Lambda_matrix)

# Check if Lambda_matrix is symmetrical
cat("\nThe matrix is symmetrical ? -", isSymmetric(Lambda_matrix), "\n") 

# Get vector of h2_intercept for each phenotype
h2_intercept = genetic_correlation %>%
  select(ID2, 'h2 intercept') %>%
  slice(1:size_matrix)

h2_intercept = h2_intercept[c(size_matrix, 1:(size_matrix-1)),]

# Update Diagonal of Lambda_matrix
diag(Lambda_matrix) = h2_intercept$`h2 intercept`

# LAMBDA ZERO MATRIX COMPLETE

### Completing the Lambda matrix with the sample size of each phenotype
# Create a vector of sample size
# Check concordance between genetic_correlation and output_loglik
same_column = output_Loglik$Phenotype == rownames(Lambda_matrix)
if (all(same_column) == FALSE) {
  stop("Rownames between genetic correlation and Log-likelihood output are not identical")
}

# retrieve sample size and take the squared root of it
sqrt_sample_size = 1/sqrt(output_Loglik$Sample_size)

# Complete Lambda matrix 
Lambda_matrix = diag(sqrt_sample_size) %*% Lambda_matrix %*% diag(sqrt_sample_size)

## LAMBDA MATRIX IS COMPLETE

###################################################
# Get squared root of h2 intercept
h2_intercept = h2_intercept %>%
  mutate(root_h2_intercept = h2_intercept$`h2 intercept`^(1/2))


###########################

# PRODUCE CORRELATION PLOT OF ALL PHENOTYPES

###########################

# # Correlation plot
# corrplot(matrix_correlation, method = "color", type = "full", order = "hclust",
#          tl.col="black", tl.srt=45, tl.cex = 0.6, tl.offset = 0.5,
#          p.mat = significance_matrix, sig.level = 0.05, insig = "blank"
#          #col = brewer.pal(n = 8, name = "PuOr")
# )

###########################

# CLUSTERING OF THE DATA

###########################
# Here the goal is to separate traits into semi-correlated traits
# And then the calculation will be performed within those cluster of semi-correlated traits
# First we observed how much cluster is appropriate
# By looking at the correlation plot we can already observe that 3 clusters are well represented


# require(vegan)
# require(FactoMineR)
# require(gclus)
# library(cluster)
# library(fpc)
# library(ClusterR)

###### K-Mean clustering, determination of optimal number of cluster
# Elbow method, 
# k.max <- 20
# 
# # we will create a vector of the total within sum of squars, in order to visualize it 
# wss <- sapply(1:k.max, function(k){kmeans(matrix_correlation, k, 
#                                           nstart=50,iter.max = 1000 )$tot.withinss})

# We can use the built in function fviz_nbclust to display "Elbow": 
# fviz_nbclust(matrix_correlation, FUN = kmeans,method = "wss" ,nstart = 50, 
#              linecolor = "red") +
#   theme(axis.text.x = element_text(size = 15),title = element_text(size = 15, color = "black"))
# plot(1:k.max, wss, type="b", pch = 19,  xlab="Number of clusters K", ylab="Total within-clusters sum of squares")
# 
# # Silhouette method
# opt.k.sil <- Optimal_Clusters_KMeans(matrix_correlation, max_clusters=20, plot_clusters = TRUE, criterion = "silhouette")
# fviz_nbclust(matrix_correlation ,kmeans , method = "silhouette",
#              linecolor = "red") +
#   theme(axis.text.x = element_text(size = 15),title = element_text(size = 15, color = "black"))
# 
# # Calinski-Harabasz method
# for (i in 2:k.max) {
#   km = kmeans(matrix_correlation, centers = i, nstart = 50)
#   value = round(calinhara(matrix_correlation, km$cluster), digits = 1)
#   cat("Number of cluster:", i, " value of Calinski_Harabasz : ", value, "\n")
# }

# The optimal number of cluster seems to be 2, 7 or 8

###########################

# CLUSTERING OF THE DATA INTO 3 LISTs

###########################

number_center = 3
km.result = kmeans(matrix_correlation, centers = number_center, nstart = 10)


# Plot of the clusters 
# fviz_cluster(km.result, data = matrix_correlation, geom = c("point"), ellipse.type = "norm",
#              ggtheme = theme_bw()) +
#   theme(axis.text.x = element_text(size = 15),title = element_text(size = 15, color = "black"))

# Create subset of cluster
cluster_list = list()

# Formation of a list of cluster containing all the subset cluster
for (i in 1:number_center) {
  x = which(km.result$cluster == i)
  cluster_list[[i]] = matrix_correlation[x,x]
  names(cluster_list)[i] = paste("cluster",i, sep = "_")
}

# Now the goal is to analysis the correlation between all these traits within each cluster
# Traits with a too strong correlation should be discarded, here the limit is 0.7
# Might also remove traits without any correlation (pruning)
treshold = 0.7
for (i in 1:length(cluster_list)) {
  cat("\nIterating the cluster number:", i, "\n")
  x = which((cluster_list[[i]] > treshold | cluster_list[[i]] < -treshold  ) & cluster_list[[i]] != 1, arr.ind = TRUE)
  
  while (nrow(x) != 0) {
    cluster_list[[i]] = cluster_list[[i]][-x[1],-x[1]]
    cat("\nSelecting out:", rownames(x)[1], "\n")
    x = which((cluster_list[[i]] > treshold | cluster_list[[i]] < -treshold ) & cluster_list[[i]] != 1, arr.ind = TRUE)
  }
}

# Remove traits with a mean correlation below a certain threshold
threshold = 0.05
for (i in 1:length(cluster_list)) {
  cat("\nIterating the cluster number:", i, "\n")
  x = NULL
  for (j in 1:nrow(cluster_list[[i]])) {
    if ((sum(cluster_list[[i]][j,])-1)/nrow(cluster_list[[i]]) < threshold) {
      x[j] = TRUE
      cat("\nSelecting out:", rownames(cluster_list[[i]])[j], "\n")
    } else {
      x[j] = FALSE
    }
  }
  cluster_list[[i]] = cluster_list[[i]][!x,!x]
}

# Now my clusters contains semi-correlated traits !

# Save name of discarded traits 
# The goal is to redo a clustering without these traits and check for difference
list_conserved_traits = list()
for (i in 1:length(cluster_list)) {
  list_conserved_traits[[i]] = rownames(cluster_list[[i]])
}
conserved_traits = unlist(list_conserved_traits) # !! change the order specifically for 3 cluster !! 
#conserved_traits
discarded_traits = rownames(matrix_correlation)[!(rownames(matrix_correlation) %in% conserved_traits)]
#discarded_traits

# Eliminate these discarded traits from the correlation matrix and check the clustering
x = which(rownames(matrix_correlation) %in% conserved_traits)
matrix_conserved_traits = (matrix_correlation[x,x])
significance_matrix_conserved_traits = significance_matrix[x,x]

# conserve order from the clustering method to get a nice correlation plot
order_string_sequence = NULL
for (i in 1:length(conserved_traits)) {
  order_string_sequence[i] = which(conserved_traits[i] == rownames(matrix_conserved_traits))
}
matrix_conserved_traits = matrix_conserved_traits[order_string_sequence, order_string_sequence]
significance_matrix_conserved_traits = significance_matrix_conserved_traits[order_string_sequence, order_string_sequence]

# Check the order of the matrix of conserved element
rownames(matrix_conserved_traits) == conserved_traits
#rownames(significance_matrix_conserved_traits) == conserved_traits

# Correlation plot of the consertved element
# corrplot(matrix_conserved_traits, method = "color", type = "full", 
#          tl.col="black", tl.srt=45, tl.cex = 0.6, tl.offset = 0.5,
#          p.mat = significance_matrix_conserved_traits  , sig.level = 0.05, insig = "blank"
#          #col = brewer.pal(n = 8, name = "PuOr")
# )

## Redo clustering 
##### K-mean clustering
# Does it make sense to do such graph ??
km.result = kmeans(matrix_conserved_traits, centers = number_center, nstart = 10)
km.result

# fviz_cluster(km.result, data = matrix_conserved_traits, geom = c("point"), ellipse.type = "norm",
#              ggtheme = theme_bw()) + 
#   theme(axis.text.x = element_text(size = 15),title = element_text(size = 15, color = "black"))


cat("First cluster contains: \n", rownames(cluster_list[[1]]), "\n")
cat("Second cluster contains: \n",rownames(cluster_list[[2]]), "\n")
cat("Third cluster contains: \n", rownames(cluster_list[[3]]), "\n")


####################################

# PREDICTOR/SNP AND CORRESPONDING LD SCORE REGRESSION

####################################

# Open file containing all SNP and corresponding LD score regression
LD_score = read.table(input3, 
                      sep = ",",
                      stringsAsFactors = FALSE, header = TRUE)

head(LD_score, n = 10)

# Use only a portion of LD_score to be more efficient
n = seq(1,nrow(LD_score), 80)
LD_score = LD_score[n,]

head(LD_score)
# Creating Omega matrix for one cluster
# Get index of conserved traits

# select cluster with 1558 in it
# Only necessary for my application..

is_cluster1 = any(rownames(cluster_list$cluster_1) == "1558")
is_cluster2 = any(rownames(cluster_list$cluster_2) == "1558")
is_cluster3 = any(rownames(cluster_list$cluster_3) == "1558")

selected_cluster = NULL
if (is_cluster1) {
  selected_cluster = 1
} else if (is_cluster2) {
  selected_cluster = 2
} else if (is_cluster3){
  selected_cluster = 3
}
cat("\nSelected cluster number: ",selected_cluster,"\n")

correlation_cluster_list = cluster_list
index_cluster = which(rownames(matrix_correlation) %in% rownames(correlation_cluster_list[[selected_cluster]]), arr.ind = TRUE)


# Creation of a tau matrix
tau = output_Loglik[index_cluster, "Tau"]
Phenotypes_cluster_name = output_Loglik[index_cluster, "Phenotype"]
names(tau) = Phenotypes_cluster_name
tau_matrix = sqrt(diag(tau) %*% matrix(1,length(tau), length(tau)) %*% diag(tau)) # squared root
colnames(tau_matrix) = Phenotypes_cluster_name
rownames(tau_matrix) = Phenotypes_cluster_name

# creation of alpha matrix
alpha = output_Loglik[index_cluster, "Alpha"]
names(alpha) = Phenotypes_cluster_name
alpha_matrix = (diag(alpha) %*% matrix(1,length(alpha), length(alpha)) + t(diag(alpha) %*% matrix(1,length(alpha), length(alpha))))/2
colnames(alpha_matrix) = Phenotypes_cluster_name
rownames(alpha_matrix) = Phenotypes_cluster_name

# Creation of matrix of correlation
correlation_matrix = correlation_cluster_list[[selected_cluster]]

# LD score function giving a vector 
# LD_score_function = function(SNP_number,alp){
#   value = LD_score$LD_Intercept[SNP_number] + LD_score$LD_linear[SNP_number] * alp + LD_score$LD_quadratic[SNP_number] * alp^2 + LD_score$LD_cubic[SNP_number] * alp^3 + LD_score$LD_quartic[SNP_number] * alp^4
#   return(value)
# }
# 
# ##### Create a matrix of LD score for each SNP where the value are stored in vector
# 
# # LD score function giving a vector 
# LD_score_function = function(alp){
#   value = LD_score$LD_Intercept + LD_score$LD_linear * alp + LD_score$LD_quadratic * alp^2 + LD_score$LD_cubic * alp^3 + LD_score$LD_quartic * alp^4
#   return(value)
# }

# declaring lower triangle matrix so that the memory usage is smaller
lower_alpha_matrix = alpha_matrix[lower.tri(alpha_matrix, diag = TRUE)]
lower_tau_matrix = tau_matrix[lower.tri(tau_matrix, diag = TRUE)]
lower_correlation_matrix = correlation_matrix[lower.tri(correlation_matrix, diag = TRUE)]
lambda_matrix_cluster = Lambda_matrix[index_cluster,index_cluster]
lower_lambda_matrix = lambda_matrix_cluster[lower.tri(lambda_matrix_cluster, diag = TRUE)]

####################################

# COMPUTE OMEGA AND LAMBDA MATRIX

####################################

# column vector containing alpha to the corresponding power
# will be used for the calculation of Li_matrix later
li_matrix_score_2 = matrix(NA,ncol = length(lower_alpha_matrix), nrow = 5)
for (i in 1:length(lower_alpha_matrix)) {
  x = c(1,lower_alpha_matrix[i], lower_alpha_matrix[i]^2, lower_alpha_matrix[i]^3, lower_alpha_matrix[i]^4)
  li_matrix_score_2[,i] = x
}

# define matrix of LD coefficient without predictor name just for calculation
temp_LD_score = as.matrix(LD_score[,-c(1,2)])

head(li_matrix_score_2)
length(li_matrix_score_2)
head(temp_LD_score)
length(temp_LD_score)

# Calculate LD Score by multiplication of coefficient and matrix of alpha to 4 degree
Omega_LD_score = temp_LD_score %*% li_matrix_score_2

# Applying multiplication of correlation matrix, tau matrix to get a full Omega matrix
Omega = Omega_LD_score %*% diag(lower_correlation_matrix) %*% diag(lower_tau_matrix)

# Adding Lambda matrix to Omega 
Omega_Lambda = Omega + rep(lower_lambda_matrix, each = nrow(Omega))



####################################

# COMPUTE PDF FUNCTION 

####################################


# package
library(MASS)

# Compute Logarithm of PDF function
Log_PDF = function(B, Sig){
  k = nrow(Sig)
  x = -(k/2)*log(2*pi) - 1/2*log(det(Sig)) - 1/2*(t(B) %*% ginv(Sig) %*% B) #solve(Sig)
  return(x)
}

# Function for producing symmetrical matrix
makeSymm <- function(m) {
  m[upper.tri(m)] = t(m)[upper.tri(m)]
  return(m)
}

# We need the Beta value for each phenotype
Beta_complete = read.table(input1, 
                           sep = ",",
                           stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

Beta_complete = Beta_complete[n,]
head(Beta_complete)

# Compare the predictors names between LD_score and Z_score
Same_predictor = LD_score$Predictor == Beta_complete$Predictor
cat("\nSame predictor between LD_score and Beta table? ", all(Same_predictor), "\n")

# Select the column of Beta that are corresponding to the selected traits
Beta = as.matrix(Beta_complete[,Phenotypes_cluster_name])
head(Beta)

# Produce Sum of log PDF function
Sum_Log = 0
Sigma_matrix = matrix(NA, nrow = nrow(tau_matrix), ncol = nrow(tau_matrix))
for (i in 1:nrow(LD_score)) {
  Sigma_matrix[lower.tri(Sigma_matrix, diag = TRUE)] = Omega_Lambda[i,]
  
  Sigma_matrix = makeSymm(Sigma_matrix)
  Beta_i = Beta[i,]
  Sum_Log = Sum_Log + Log_PDF(Beta_i, Sigma_matrix)
}

#######################################################################

Big_PDF_function = function(alpha_tau, Lambda_par, correlation, LD, Beta_score){
  
  alpha_par = alpha_tau[1:(length(alpha_tau)/2)]
  tau_par = alpha_tau[(length(alpha_tau)/2+1):length(alpha_tau)]
  
  alpha_complete = (diag(alpha_par) %*% matrix(1,length(alpha_par), length(alpha_par)) + t(diag(alpha_par) %*% matrix(1,length(alpha_par), length(alpha_par))))/2
  lower_alpha = alpha_complete[lower.tri(alpha_complete, diag = TRUE)]
  
  tau_complete = sqrt(diag(tau_par) %*% matrix(1,length(tau_par), length(tau_par)) %*% diag(tau_par))
  lower_tau = abs(tau_complete[lower.tri(tau_complete, diag = TRUE)])
  
  lower_lambda = Lambda_par
  lower_correlation = correlation
  LD_score_PDF = LD
  
  # determine value of weighted LD score - to optimize
  li_matrix_score = matrix(NA,ncol = length(lower_alpha), nrow = 5)
  for (i in 1:length(lower_alpha)) {
    x = c(1,lower_alpha[i], lower_alpha[i]^2, lower_alpha[i]^3, lower_alpha[i]^4)
    li_matrix_score[,i] = x
  }
  
  # define matrix of LD coefficient without predictor name just for calculation
  temp_LD_score_PDF = as.matrix(LD_score_PDF[,-c(1,2)]) 
  
  # Calculate LD Score by multiplication of coefficient and matrix of alpha to 4 degree
  Omega_LD_score_PDF = temp_LD_score_PDF %*% li_matrix_score
  
  # Applying multiplication of correlation matrix, tau matrix to get a full Omega matrix
  Omega_PDF = Omega_LD_score_PDF %*% diag(lower_correlation) %*% diag(abs(lower_tau))
  
  # Adding Lambda matrix to Omega 
  Omega_Lambda_PDF = Omega_PDF + rep(lower_lambda, each = nrow(Omega_PDF))
  Omega_Lambda_PDF[is.na(Omega_Lambda_PDF)] = 0
  
  #cat("\nNA in OMEGA LAMBDA ?  ",any(is.na(Omega_Lambda_PDF)))
  #cat("\nInfinite value ?  ", any(is.infinite(Omega_Lambda_PDF)), "\n")
  
  size = ncol(Beta_score)
  
  Sum_Log_PDF = 0
  Sigma_matrix_PDF = matrix(NA, nrow = size, ncol = size)
  for (i in 1:nrow(LD_score_PDF)) {
    Sigma_matrix_PDF[lower.tri(Sigma_matrix_PDF, diag = TRUE)] = Omega_Lambda_PDF[i,]
    
    Sigma_matrix_PDF = makeSymm(Sigma_matrix_PDF)
    Beta_i = Beta[i,]
    Sum_Log_PDF = Sum_Log_PDF + Log_PDF(Beta_i, Sigma_matrix_PDF)
  }
  
  cat(alpha_par[1:3], "** **", tau_par[1:3], "** **", Sum_Log_PDF, "\n")
  
  return(Sum_Log_PDF)
}

parameter_together = c(alpha,tau)

Big_PDF_function(parameter_together, lower_lambda_matrix, lower_correlation_matrix, 
                 LD_score, Beta)


#######################################################################
# OPTIMIZATION ALL TOGETHER

print(alpha)
print(tau)

start_time = Sys.time()
print("***************************************************** OPTIMIZATION OF TAU AND ALPHA")

optim_complete = optim(par = parameter_together, fn = Big_PDF_function, 
                  Lambda_par = lower_lambda_matrix, 
                  correlation = lower_correlation_matrix, 
                  LD = LD_score, 
                  Beta_score = Beta, 
                  method = input4,
                  control = list(fnscale = -1, maxit = 200000, reltol = 1e-4)) 
end_time = Sys.time()
print(end_time - start_time)


cat("\n********************************\n", "END OF OPTIMIZATION","\n********************************\n" )

print(input4)

cat("Optim value")
optim_complete$value

cat("\n alpha optimized")
length_optim = length(optim_complete$par)
optim_complete$par[1:(length_optim/2)]
cat("\n tau optimized")
optim_complete$par[(length_optim/2+1):length_optim]

# The optimized parameters need to be extracted from the log files at this point






