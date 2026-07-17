

COMPUTE_BGE_FC <- function(T_m_LIST,m_LIST,incidence,a,v,T_0){
  
  n = nrow(incidence)

  P_local_num = rep(0,n) 
  P_local_den = rep(0,n) 
  
  

  for (i_data in 1:length(T_m_LIST)){
    
    T_m_i = T_m_LIST[[i_data]]
    m_i   = m_LIST[[i_data]]
 
    if(m_i>0){
    
      for (j in 1:n)  {
        n_nodes <- which(incidence[,j]==1)         # parents of j
        P_local_num[j]   <-  P_local_num[j] + (-(length(n_nodes)+1)*m_i/2)*log(2*pi) + ((length(n_nodes)+1)/2)*log(v/(v+m_i)) + c_function((length(n_nodes)+1),a-n+length(n_nodes)+1)  - c_function((length(n_nodes)+1),a+m_i-n+length(n_nodes)+1)        + ((a-n+length(n_nodes)+1)/2)*log(det(as.matrix(T_0[sort(c(n_nodes,j)),sort(c(n_nodes,j))])))   + (-(a-n+length(n_nodes)+1+m_i)/2) *log(det(as.matrix(T_m_i[sort(c(n_nodes,j)),sort(c(n_nodes,j))])))
        if(sum(incidence[,j])>0){          # if j has at least one parent
          P_local_den[j] <-  P_local_den[j] + (-(length(n_nodes))*m_i/2)*log(2*pi) + (length(n_nodes)/2)*log(v/(v+m_i))     + c_function(length(n_nodes),    a-n+length(n_nodes)  )  - c_function(length(n_nodes),    a+m_i-n+length(n_nodes)  )        + ((a-n+length(n_nodes)  )/2)*log(det(as.matrix(T_0[n_nodes,n_nodes])))                         + (-(a-n+length(n_nodes)  +m_i)/2) *log(det(as.matrix(T_m_i[n_nodes,n_nodes])))
        }
      }
    }
  }

  bge_old = (sum(P_local_num))-(sum(P_local_den))
  
  return(list(bge_old,P_local_num,P_local_den))
  
}