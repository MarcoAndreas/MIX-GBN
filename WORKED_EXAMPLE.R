
###############################################################################
# Installation and setup
###############################################################################

# 1. Save all R files in a folder named, for example, mGBN. 
#
# 2. Set the working directory to that folder. 
#    Replace X: with the appropriate drive letter on your computer. 

setwd("X:/mGBN") 

# Install the required CRAN packages (only needed once)

install.packages(c("pracma", "MASS", "MCMCpack"), repos = "https://cloud.r-project.org") 

# Load the packages 

library(pracma) 
library(MASS) 
library(MCMCpack)

# Only a few basic functions from these packages are used.

######################################################################
# Load user-defined functions for Bayesian networks and structure MCMC
######################################################################

source("generic_functions.R")

#############################################
# Functions for the homogeneous GBN model (H)
#############################################

source("SAMPLE_PARAMETERS_H.R")

source("COMPUTE_BGE_H.R")

source("COMPUTE_BGE_LOCAL_H.R")

source("MCMC_H.R")

source("PREDICTIVE_H.R")

###########################################################
# Functions for the full-covariance mixture GBN model (FC)
###########################################################

source("COMPUTE_BGE_FC.R")

source("COMPUTE_BGE_LOCAL_FC.R")

source("MCMC_FC.R")

source("GENERATE_PARA_LIST_FC.R")

source("PREDICTIVE_FC.R")

source("SAMPLE_PARAMETERS_FC.R")

##########################################################
# Functions for the tied-covariance mixture GBN model (TC)
##########################################################

source("COMPUTE_BGE_TC.R")

source("COMPUTE_BGE_LOCAL_TC.R")

source("MCMC_TC.R")

source("GENERATE_PARA_LIST_TC.R")

source("PREDICTIVE_TC.R")

source("SAMPLE_PARAMETERS_TC.R")

################################
# Example: Fisher's Iris dataset
################################

# We demonstrate the usage of the implemented methods 
# using Fisher's Iris dataset, which is included in the 
# standard R installation.

# Prepare and preprocess the Iris data.

# Extract the data matrix.
data = as.matrix(iris[, 1:4])

# Number of variables (nodes) and observations.
n = ncol(data)   # n = 4   
m = nrow(data)   # m = 150

# The n=4 variables are
# X1: Sepal.Length
# X2: Sepal.Width 
# X3: Petal.Length 
# X4: Petal.Width 

# Standardize each variable to have mean 0 and variance 1.
for (i in 1:n) {
  data[, i] = (data[, i] - mean(data[, i])) / sd(data[, i])
}

###############################################################################
# MCMC simulation settings
###############################################################################

# For illustration, we use very short MCMC runs.
# For a more thorough analysis, uncomment the settings below.

iterations = 5000
step_save  =    5

# Recommended settings for the analyses in the paper:
# Warning: Running the analyses may take about one hour.
# iterations = 500000
# step_save  =    500

# Interpretation:
# iterations = total number of MCMC iterations
# step_save  = save every step_save-th MCMC sample (thinning interval)

######################################################################
# Initial values
######################################################################

# Initial directed acyclic graph (DAG).

DAG = matrix(0, n, n)

# Initial allocation vector 
# (with three randomly initialized mixture components).

VEC_DIS = sample(1:10, m, replace = TRUE)

########################################################################
# Run the MCMC algorithm for each of the three models
########################################################################

# Homogeneous GBN model (H).
out_H <- MCMC_H(data, DAG, iterations, step_save)

# Full-covariance mixture GBN model (FC).
out_FC <- MCMC_FC(data, VEC_DIS, DAG, iterations, step_save)

# Tied-covariance mixture GBN model (TC).
out_TC <- MCMC_TC(data, VEC_DIS, DAG, iterations, step_save)

###############################################################################
# Inspect the sampled DAGs
###############################################################################

out_H[[1]]
out_FC[[1]]
out_TC[[1]]

###############################################################################
# Marginal edge posterior probabilities
###############################################################################

# Discard the first 50% of the sampled DAGs as burn-in.

n_burn = floor(0.5 * iterations / step_save)

# Extract the CPDAGs corresponding to the sampled DAGs after burn-in.

CPDAG_LIST_H  = cpdag_list(out_H[[1]],  n_burn + 1)
CPDAG_LIST_FC = cpdag_list(out_FC[[1]], n_burn + 1)
CPDAG_LIST_TC = cpdag_list(out_TC[[1]], n_burn + 1)

# Posterior edge probabilities.

CPDAG_H  = CPDAG_LIST_H[[3]]
CPDAG_FC = CPDAG_LIST_FC[[3]]
CPDAG_TC = CPDAG_LIST_TC[[3]]

CPDAG_H
CPDAG_FC
CPDAG_TC

# Example:

CPDAG_TC[1,3] 

# is the marginal edge posterior probability of the directed edge
# X1 (Sepal.Length) -> X3 (Petal.Length)
# (inferred with the tied-covariance mixture GBN model (TC))



#########################################################################
# Compute the co-allocation probabilities for the TC results
#########################################################################

matrix_heat = co_allocation_matrix(out_TC)
# The first 50% of the sampled allocation vectors are discarded as burn-in.

# Display the co-allocation probabilities as a heatmap
heatmap(matrix_heat,Rowv = NA,Colv = NA, scale = "none", col = gray(seq(1, 0, length = 256)))

######################################################################## 
# Compute log geometric mean posterior predictive probabilities 
########################################################################

# Choose 100 observations for training. 

m_train = 100 

# Randomly split the observations into training and validation sets. 
ind = sample(1:m, m, replace = FALSE) 

ind_train = ind[1:m_train] 
ind_valid = ind[(m_train + 1):m] 

# Generate the training and validation datasets. 

data_train = data[ind_train, , drop = FALSE] 
data_valid = data[ind_valid, , drop = FALSE] 

m_valid = nrow(data_valid)


# Randomly initialize an allocation vector of length m_train. 
VEC_DIS = sample(1:3, m_train, replace = TRUE) 


# For each of the three models, run an MCMC simulation on the training data
# and compute the log geometric mean posterior predictive probability for the
# validation data.

#############################################################################
# Homogeneous GBN model (H)
#############################################################################
 
out_H2 = MCMC_H(data_train, DAG, iterations, step_save)
 
MATRIX = PREDICTIVE_H(out_H2, data_valid, iterations, step_save, n)
 
log_scores = rowSums(MATRIX)
 
# Log of geometric mean posterior predictive probability
# (numerically stabilized computation).
 
pred_H = (log(mean(exp(log_scores - max(log_scores)))) + max(log_scores)) / m_valid
 
#############################################################################
# Full-covariance mixture GBN model (FC)
#############################################################################
 
out_FC2 = MCMC_FC(data_train, VEC_DIS, DAG, iterations, step_save)
 
MATRIX = PREDICTIVE_FC(out_FC2, data_train, data_valid)
 
log_scores = rowSums(MATRIX)
 
# Log of geometric mean posterior predictive probability
# (numerically stabilized computation).

pred_FC = (log(mean(exp(log_scores - max(log_scores)))) + max(log_scores)) / m_valid

#############################################################################
# Tied-covariance mixture GBN model (TC)

#############################################################################
 
out_TC2 = MCMC_TC(data_train, VEC_DIS, DAG, iterations, step_save)
 
MATRIX = PREDICTIVE_TC(out_TC2, data_train, data_valid)
 
log_scores = rowSums(MATRIX)
 
# Log of geometric mean posterior predictive probability
# (numerically stabilized computation).
 
pred_TC = (log(mean(exp(log_scores - max(log_scores)))) + max(log_scores)) / m_valid
 
#########################################################################
# Display the log geometric mean posterior predictive probabilities.
#########################################################################
 
vec = c(H = pred_H, FC = pred_FC, TC = pred_TC)

vec

#####################################################################
# End of worked example
#####################################################################
 












