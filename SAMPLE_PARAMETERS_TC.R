

SAMPLE_PARAMETERS_TC <- function(T_m_LIST, m_LIST, MEANS_LIST, a, v, mu, incidence, m){
  
  n = nrow(incidence)
  
  COV_LIST = NULL
  MU_LIST  = NULL
  
  W = rWishart(1,a+m,solve(T_m_LIST[[1]]))
  W = W[,,1]
  W = (W + t(W))/2
  
  R       = chol(W) 
  R_inv   = solve(R,diag(n))   
  SIGMA   = R_inv %*% t(R_inv)
  
  SIGMA_NEW = EXTRACT_PRECISION_OF_DAG_NEW(SIGMA,incidence)
  

  for (i_data in 1:length(m_LIST)){
    
    m_i    = m_LIST[[i_data]]
    MEAN_i = MEANS_LIST[[i_data]]
    
    if(m_i>0){MU_m_i = (v*mu  + m_i*MEAN_i)/(m_i+v)}
         else{MU_m_i = mu} 
    
    MU_NEW_i =    mvrnorm(1,MU_m_i,SIGMA/(m_i+v) ) 
    
    MU_LIST[[i_data]]  = MU_NEW_i
    
  }

  COV_LIST[[1]] = SIGMA_NEW
  return(list(MU_LIST,COV_LIST))
  
}