

SAMPLE_PARAMETERS_H <- function(T_m,m,data,a,v,mu,incidence){
  
  COV_LIST = NULL
  MU_LIST = NULL
  
  n = nrow(incidence)
  
    W = rWishart(1,a+m,solve(T_m))
    W = W[,,1]
    W = (W + t(W))/2
    
    R       = chol(W) 
    R_inv   = solve(R,diag(n))   
    SIGMA   = R_inv %*% t(R_inv)
    
    
    SIGMA_NEW = EXTRACT_PRECISION_OF_DAG_NEW(SIGMA,incidence)
    
    MU_m = (v*mu  + m*rowMeans(data))/(m+v)
         
    MU_NEW = mvrnorm(1,MU_m,SIGMA/(m+v)) 
    
    COV_LIST[[1]] = SIGMA_NEW
    MU_LIST[[1]]  = MU_NEW
    
  return(list(MU_LIST,COV_LIST))
  
}