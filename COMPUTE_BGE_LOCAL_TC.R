

COMPUTE_BGE_LOCAL_TC <- function(operation,child,parent,T_m_LIST,m_LIST,m,P_local_num,P_local_den,incidence_new,a,v,T_0){
  
  n = nrow(incidence_new)
  
  T_m = T_m_LIST[[1]]
  
  m_total = m
  
  P_local_num_new <- P_local_num
  P_local_den_new <- P_local_den
  
  n_nodes_new <- which(incidence_new[,child]==1)

  
  P_local_num_new[child] <- (-(length(n_nodes_new)+1)*m_total/2)*log(2*pi) + c_function((length(n_nodes_new)+1),a-n+length(n_nodes_new)+1) -c_function((length(n_nodes_new)+1),a+m_total-n+length(n_nodes_new)+1)     + ((a-n+length(n_nodes_new)+1)/2)  *log(det(as.matrix(T_0[sort(c(n_nodes_new,child)),sort(c(n_nodes_new,child))])))+ (-(a-n+length(n_nodes_new)+1+m_total)/2)  *log(det(as.matrix(T_m[sort(c(n_nodes_new,child)),sort(c(n_nodes_new,child))])))
  
  if(sum(incidence_new[,child])>0){           # if child at least one parent
    P_local_den_new[child] <- (-(length(n_nodes_new))*m_total/2)*log(2*pi)+ c_function(length(n_nodes_new)    ,a-n+length(n_nodes_new)  ) - c_function(length(n_nodes_new),   a+m_total-n+length(n_nodes_new)  )      + ((a-n+length(n_nodes_new))  /2)  *log(det(as.matrix(T_0[n_nodes_new,n_nodes_new])))                              + (-(a-n+length(n_nodes_new)+m_total)/2)  *log(det(as.matrix(T_m[n_nodes_new,n_nodes_new])))
  }
  else{P_local_den_new[child] <- 0}
  
  
  for (i_data in 1:length(m_LIST)){
    
    m_i   = m_LIST[[i_data]]
    
    P_local_num_new[child] <- P_local_num_new[child] + ((length(n_nodes_new)+1)/2)*log(v/(v+m_i)) 
    
    if(sum(incidence_new[,child])>0){           # if child at least one parent
      P_local_den_new[child] <-   P_local_den_new[child] + (length(n_nodes_new)/2)*log(v/(v+m_i))     
    }
  }

  
  
    if (operation==1){                          # if single edge operation was an edge reversal
      
      n_nodesP <- which(incidence_new[,parent]==1)
      
      P_local_num_new[parent] <- (-(length(n_nodesP)+1)*m_total/2)*log(2*pi)      + c_function((length(n_nodesP)+1),a-n+length(n_nodesP)+1) - c_function((length(n_nodesP)+1),a+m_total-n+length(n_nodesP)+1)     + ((a-n+length(n_nodesP)+1)/2)*log(det(as.matrix(T_0[sort(c(n_nodesP,parent)),sort(c(n_nodesP,parent))])))+ (-(a+m_total-n+length(n_nodesP)+1)/2)  *log(det(as.matrix(T_m[sort(c(n_nodesP,parent)),sort(c(n_nodesP,parent))])))
      if(sum(incidence_new[,parent])>0){          # if parent at least one parent
        P_local_den_new[parent] <- (-(length(n_nodesP))*m_total/2)*log(2*pi)      + c_function(length(n_nodesP)    ,a-n+length(n_nodesP)  ) - c_function(length(n_nodesP)    ,a+m_total-n+length(n_nodesP)  )     + ((a-n+length(n_nodesP)) /2)*log(det(as.matrix(T_0[n_nodesP,n_nodesP])))                                + (-(a+m_total-n+length(n_nodesP))  /2)  *log(det(as.matrix(T_m[n_nodesP,n_nodesP])))
      }
      else{P_local_den_new[parent] <- 0}
      
      for (i_data in 1:length(m_LIST)){
        
        m_i                     = m_LIST[[i_data]]
        P_local_num_new[parent] = P_local_num_new[parent]  + ((length(n_nodesP)+1)/2)*log(v/(v+m_i)) 
      if(sum(incidence_new[,parent])>0){          # if parent at least one parent
        P_local_den_new[parent] =  P_local_den_new[parent] + (length(n_nodesP)/2)*log(v/(v+m_i))   
      }
     
  
      }
    }
    

    bge_new <- (sum(P_local_num_new))-(sum(P_local_den_new))
    

  
  
  return(list(bge_new,P_local_num_new,P_local_den_new))
  
  
}