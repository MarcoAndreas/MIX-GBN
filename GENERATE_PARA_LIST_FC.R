

GENERATE_PARA_LIST_FC <- function(indicis,data_con,T_0,v,mu,m){
  
n = ncol(data_con)
  
T_m_LIST   = NULL
m_LIST     = NULL
MEANS_LIST = NULL

K = max(indicis)

for (group_index in 1:K){
  
  ind_i = which(indicis==group_index)
  
  data_i             = data_con[ind_i,, drop = FALSE]
  
  if(isempty(data_i)==FALSE){
    m_i = nrow(data_i)
    MEANS_LIST[[group_index]] = colMeans(data_i)
    
    if(nrow(data_i)>1){ 
      T_m_LIST[[group_index]] = T_0 + (m_i-1)* cov(data_i) + ((v*m_i)/(v+m_i))* (mu - MEANS_LIST[[group_index]])%*%t(mu - MEANS_LIST[[group_index]])}
    
    
    else {T_m_LIST[[group_index]] = T_0  + ((v*m_i)/(v+m_i))* (mu - MEANS_LIST[[group_index]])%*%t(mu - MEANS_LIST[[group_index]])}
  }   
  
  
  
  
  else{m_i = 0
  MEANS_LIST[[group_index]] = rep(0,n)
  T_m_LIST[[group_index]]  =  T_0
  }
  
  m_LIST[[group_index]]     = m_i
 
  
}



return(list(MEANS_LIST,T_m_LIST,m_LIST))

}
