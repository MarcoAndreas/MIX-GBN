
PREDICTIVE_FC <- function(mix_hybrid_new,data_con,data_con2){

n  <- ncol(data_con)    # number of continuous nodes
m  <- nrow(data_con)    # number of observations
m2 <- nrow(data_con2)   # number of validation observations

v=1
  
mu = numeric(n)
a  = n+2 

T_0 = diag(1,n,n)
T_0 = (a-n-1) * T_0 # scaling T_0
     
############################################################################

n_samples = length(mix_hybrid_new[[1]])

MATRIX = NULL

for (i_sample in round((n_samples/2)+1):n_samples){

DAG = mix_hybrid_new[[1]][[i_sample]]
VEC = mix_hybrid_new[[2]][[i_sample]]

pi_vec = rdirichlet(1,table(VEC)+1)

OUT = GENERATE_PARA_LIST_FC(VEC,data_con,T_0,v,mu)

MEANS_LIST = OUT[[1]]
T_m_LIST   = OUT[[2]]
m_LIST     = OUT[[3]]

out = SAMPLE_PARAMETERS_FC(T_m_LIST, m_LIST, MEANS_LIST, a, v, mu, DAG)

MU_LIST  = out[[1]]
COV_LIST = out[[2]]

K = max(VEC)

LOG_SCORES = rep(0,m2)

for (i_obs in 1:m2){

        vec_con2 = data_con2[i_obs,]

        vec_log_components = rep(0,K)

	for (k in 1:K){

  		MU    = MU_LIST[[k]]
  		SIGMA = COV_LIST[[k]]
      
    
      		R     = chol(SIGMA) 
      		R_inv = solve(R,diag(n))   
      		W     = R_inv %*% t(R_inv)
      
      		vec_log_components[k] = log(pi_vec[k]) -n/2 * log(2*pi)  +0.5 * log(det(W))  -0.5 * t(vec_con2 - MU)%*% W %*% (vec_con2 - MU)
      	}

     const = max(vec_log_components)

     LOG_SCORES[i_obs] = log(sum(exp(vec_log_components-const))) +  const

     
}

MATRIX = rbind(MATRIX,LOG_SCORES)

}

return(MATRIX)

}
 




