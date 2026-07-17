
PREDICTIVE_H <- function(out_bge,data_con2,iterations,step_save,n_nodes){

m2 = nrow(data_con2)

MATRIX  = NULL

#################################################################################


for (k in ((0.5 * iterations/step_save) + 1):length(out_bge[[3]])){ 
  
  LOG_SCORES = rep(0,m2)


  for (i_obs in 1:m2){
  

    MU    = out_bge[[3]][[k]][[1]]
    SIGMA = out_bge[[4]][[k]][[1]]
    
    R     = chol(SIGMA) 
    R_inv = solve(R,diag(n_nodes))   
    W     = R_inv %*% t(R_inv)
    
    log_p =  -n_nodes/2 * log(2*pi)  +0.5 * log(det(W))  -0.5 * t(data_con2[i_obs,] - MU)%*% W %*% (data_con2[i_obs,] - MU)

    LOG_SCORES[i_obs] = log_p
    
  }
  
   MATRIX     = rbind(MATRIX,LOG_SCORES)
   
}


return(MATRIX)
}