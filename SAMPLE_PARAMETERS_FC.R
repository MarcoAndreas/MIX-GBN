

SAMPLE_PARAMETERS_FC <- function(T_m_LIST,m_LIST,MEANS_LIST,a,v,mu,incidence){
  
  n = nrow(incidence)
  
  COV_LIST = NULL
  MU_LIST  = NULL
  
  for (i_data in 1:length(T_m_LIST)){
    
    T_m_i  = T_m_LIST[[i_data]]
    m_i    = m_LIST[[i_data]]
    MEAN_i = MEANS_LIST[[i_data]]
    

    W_i = rWishart(1,a+m_i,solve(T_m_i))
    W_i = W_i[,,1]
    W_i = (W_i + t(W_i))/2
    
    R         = chol(W_i) 
    R_inv     = solve(R,diag(n))   
    SIGMA_i   = R_inv %*% t(R_inv)
    
    
    SIGMA_NEW_i = EXTRACT_PRECISION_OF_DAG_NEW(SIGMA_i,incidence)
    
    
    MU_m_i = (v*mu  + m_i*MEAN_i)/(m_i+v)
    
    # if(m_i>0){}
         #else{MU_m_i = mu} 
    
    MU_NEW_i =    mvrnorm(1,MU_m_i,SIGMA_i/(m_i+v) ) 
    
    COV_LIST[[i_data]] = SIGMA_NEW_i
    MU_LIST[[i_data]]  = MU_NEW_i
    
  }

  
  return(list(MU_LIST,COV_LIST))
  
}