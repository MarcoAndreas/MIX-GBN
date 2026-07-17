

COMPUTE_BGE_LOCAL_FC <- function(operation, child, parent, T_m_LIST, m_LIST, P_local_num, P_local_den, incidence_new, a, v, T_0){

  n = nrow(incidence_new)
  
  P_local_num_new <- P_local_num
  P_local_den_new <- P_local_den

  n_nodes_new <- which(incidence_new[,child]==1)
  
  
  P_local_num_new[child] = 0
  P_local_den_new[child] = 0
  
  
  if (operation==1){                          # if single edge operation was an edge reversal
  
    P_local_num_new[parent] = 0
    P_local_den_new[parent] = 0
  
  }
  
  
  
  
  
  
  
  for (i_data in 1:length(T_m_LIST)){
    
    T_m_i = T_m_LIST[[i_data]]
    m_i   = m_LIST[[i_data]]
    
    if(m_i>0){
    
    P_local_num_new[child] <- P_local_num_new[child] + (-(length(n_nodes_new)+1)*m_i/2)*log(2*pi) + ((length(n_nodes_new)+1)/2)*log(v/(v+m_i)) + c_function((length(n_nodes_new)+1),a-n+length(n_nodes_new)+1) -c_function((length(n_nodes_new)+1),a+m_i-n+length(n_nodes_new)+1) + ((a-n+length(n_nodes_new)+1)/2)  *log(det(as.matrix(T_0[sort(c(n_nodes_new,child)),sort(c(n_nodes_new,child))])))+ (-(a-n+length(n_nodes_new)+1+m_i)/2)  *log(det(as.matrix(T_m_i[sort(c(n_nodes_new,child)),sort(c(n_nodes_new,child))])))
    
    
    if(sum(incidence_new[,child])>0){           # if child at least one parent
      P_local_den_new[child] <-P_local_den_new[child] + (-(length(n_nodes_new))*m_i/2)*log(2*pi) + (length(n_nodes_new)/2)*log(v/(v+m_i))     + c_function(length(n_nodes_new)    ,a-n+length(n_nodes_new)  ) - c_function(length(n_nodes_new),   a+m_i-n+length(n_nodes_new)  ) + ((a-n+length(n_nodes_new))  /2)  *log(det(as.matrix(T_0[n_nodes_new,n_nodes_new])))                              + (-(a-n+length(n_nodes_new)+m_i  )/2)  *log(det(as.matrix(T_m_i[n_nodes_new,n_nodes_new])))
   
  
      
       }
   
    }
 
    
    
    if (operation==1){                          # if single edge operation was an edge reversal
      
  
      n_nodesP <- which(incidence_new[,parent]==1)
      
      if(m_i>0){
      
      P_local_num_new[parent] <- P_local_num_new[parent] + (-(length(n_nodesP)+1)*m_i/2)*log(2*pi) + ((length(n_nodesP)+1)/2)*log(v/(v+m_i)) + c_function((length(n_nodesP)+1),a-n+length(n_nodesP)+1) - c_function((length(n_nodesP)+1),a+m_i-n+length(n_nodesP)+1)     + ((a-n+length(n_nodesP)+1)/2)*log(det(as.matrix(T_0[sort(c(n_nodesP,parent)),sort(c(n_nodesP,parent))])))+ (-(a+m_i-n+length(n_nodesP)+1)/2)  *log(det(as.matrix(T_m_i[sort(c(n_nodesP,parent)),sort(c(n_nodesP,parent))])))
      if(sum(incidence_new[,parent])>0){          # if parent at least one parent
        P_local_den_new[parent] <- P_local_den_new[parent] + (-(length(n_nodesP))*m_i/2)*log(2*pi) + (length(n_nodesP)/2)*log(v/(v+m_i))     + c_function(length(n_nodesP)    ,a-n+length(n_nodesP)  ) - c_function(length(n_nodesP)    ,a+m_i-n+length(n_nodesP)  )     + ((a-n+length(n_nodesP)) /2)*log(det(as.matrix(T_0[n_nodesP,n_nodesP])))                                + (-(a+m_i-n+length(n_nodesP))  /2)  *log(det(as.matrix(T_m_i[n_nodesP,n_nodesP])))
      }
     
    }
  
    }
  
   
  }
  
  bge_new <- (sum(P_local_num_new))-(sum(P_local_den_new))
  

  
  return(list(bge_new,P_local_num_new,P_local_den_new))
  
}