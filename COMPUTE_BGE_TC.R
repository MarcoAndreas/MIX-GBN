

COMPUTE_BGE_TC <- function(T_m_LIST, m_LIST, m, incidence, a, v, T_0){
  
  n = nrow(incidence)

  P_local_num = rep(0,n) 
  P_local_den = rep(0,n) 
  
  T_m = T_m_LIST[[1]]
  
  m_total = m
  
  n = nrow(incidence)
  
  for (j in 1:n)  {
    n_nodes <- which(incidence[,j]==1)         # parents of j
    P_local_num[j] <- (-(length(n_nodes)+1)*m_total/2)*log(2*pi)    + c_function((length(n_nodes)+1),a-n+length(n_nodes)+1)  - c_function((length(n_nodes)+1),a+m_total-n+length(n_nodes)+1)        + ((a-n+length(n_nodes)+1)/2)*log(det(as.matrix(T_0[sort(c(n_nodes,j)),sort(c(n_nodes,j))])))   + (-(a-n+length(n_nodes)+1+m_total)/2) *log(det(as.matrix(T_m[sort(c(n_nodes,j)),sort(c(n_nodes,j))])))
    if(sum(incidence[,j])>0){          # if j has at least one parent
      P_local_den[j] <- (-(length(n_nodes))*m_total/2)*log(2*pi)    + c_function(length(n_nodes),    a-n+length(n_nodes)  )  - c_function(length(n_nodes),    a+m_total-n+length(n_nodes)  )        + ((a-n+length(n_nodes)  )/2)*log(det(as.matrix(T_0[n_nodes,n_nodes])))                         + (-(a-n+length(n_nodes)  +m_total)/2) *log(det(as.matrix(T_m[n_nodes,n_nodes])))
    }
    else{                              # if j has no parents
      P_local_den[j] <- 0
    }
  }



  for (i_data in 1:length(m_LIST)){
    
    m_i   = m_LIST[[i_data]]
    
      for (j in 1:n)  {
        n_nodes <- which(incidence[,j]==1)         # parents of j
        P_local_num[j] <-  P_local_num[j] + ((length(n_nodes)+1)/2)*log(v/(v+m_i))         
        if(sum(incidence[,j])>0){          # if j has at least one parent
          P_local_den[j] <-  P_local_den[j] + (length(n_nodes)/2)*log(v/(v+m_i))}
      }
    }
    
  
  bge_old = (sum(P_local_num))-(sum(P_local_den))

  return(list(bge_old,P_local_num,P_local_den))
  
}