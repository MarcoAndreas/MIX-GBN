

GENERATE_PARA_LIST_TC <- function(indicis,data_con,T_0,v,mu){
  
n = ncol(data_con)
  
T_m_LIST   = NULL
m_LIST     = NULL
MEANS_LIST = NULL

T_m_LIST[[1]]  = T_0

K = max(indicis)

for (group_index in 1:K){
  
  ind_i  = which(indicis==group_index) 
  data_i = data_con[ind_i,, drop = FALSE]
  m_i    = nrow(data_i)
  
  MEANS_LIST[[group_index]] = colMeans(data_i)
 
    if(nrow(data_i)>1){T_m_LIST[[1]] = T_m_LIST[[1]] + (m_i-1)* cov(data_i) + ((v*m_i)/(v+m_i))* (mu - MEANS_LIST[[group_index]])%*%t(mu - MEANS_LIST[[group_index]])}
                 else {T_m_LIST[[1]] = T_m_LIST[[1]]                        + ((v*m_i)/(v+m_i))* (mu - MEANS_LIST[[group_index]])%*%t(mu - MEANS_LIST[[group_index]])}
  
  m_LIST[[group_index]] = m_i
  
}


return(list(MEANS_LIST,T_m_LIST,m_LIST))

}

