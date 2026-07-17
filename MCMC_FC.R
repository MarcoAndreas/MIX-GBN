
MCMC_FC <- function(data_con, VEC_DIS, incidence,iterations,step_save){
 
  n  <- ncol(data_con)    # number of continuous nodes
  m  <- nrow(data_con)    # number of observations
  
  fan.in = n-1
  
  v=1
  
  mu=numeric(n)
  
  a=n+2 

  T_0 = diag(1,n,n)
  T_0 = (a-n-1) * T_0 # scaling T_0
  
  L1 <- list()    # DAGs
  L2 <- list()    # VEC_DIS
  L3 <- list()    # log score

################################################################################
### first graph
 
OUT = GENERATE_PARA_LIST_FC(VEC_DIS,data_con,T_0,v,mu)
  
MEANS_LIST = OUT[[1]]
T_m_LIST   = OUT[[2]]
m_LIST     = OUT[[3]]
  
################################################################################
  
  out = COMPUTE_BGE_FC(T_m_LIST, m_LIST, incidence, a, v, T_0)
 
  bge_old     = out[[1]]
  P_local_num = out[[2]]
  P_local_den = out[[3]]

  L1[[1]] <- incidence                    
  L2[[1]] <- VEC_DIS       
  L3[[1]] <- bge_old + compute_log_prior_allocation(VEC_DIS) + log(dpois(max(VEC_DIS), 1))
  
  # first ancestor matrix
  ancestMat <- ancestor(incidence)
  
  ####################################################################################################################################
 
  for (z in 2:((iterations/step_save)+1)){
    
    for (count in 1:step_save){
      
      out = PROPOSE_NEW_DAG(incidence,ancestMat,fan.in)
      
      incidence_new  = out[[1]]
      ancestMat_new  = out[[2]]
      p_forward      = out[[3]]
      p_backward     = out[[4]]
      operation      = out[[5]]
      child_node     = out[[6]]
      parent_node    = out[[7]]
      
   
      out = COMPUTE_BGE_LOCAL_FC(operation, child_node, parent_node, T_m_LIST, m_LIST, P_local_num, P_local_den, incidence_new, a, v, T_0)
      
      bge_new         = out[[1]]
      P_local_num_new = out[[2]]
      P_local_den_new = out[[3]]
      
      acceptance <- exp((bge_new +log(p_backward)) - (bge_old +log(p_forward)))
      rand       <- runif(1)
      
      if(acceptance > rand){
        incidence    <- incidence_new
        ancestMat    <- ancestMat_new
        bge_old      <- bge_new
        P_local_num  <- P_local_num_new
        P_local_den  <- P_local_den_new

      }
      
    

 VEC_DIS_NEW = VEC_DIS

      K = max(VEC_DIS_NEW)

      r = sample(1:K,1)

      ind_r = which(VEC_DIS_NEW==r)

      if (length(ind_r)==1){
     		K = K-1
     		ind_decrease              = which(VEC_DIS_NEW>r)
     		VEC_DIS_NEW[ind_decrease] = VEC_DIS_NEW[ind_decrease] - 1
		r_out                     = ind_r[1]    
     }  else{r_out = sample(ind_r,1)}
                   
     VEC_DIS_NEW[r_out]     = 0


    SCORES_VEC = rep(0,K+1)
    
    for (k in 1:K){
     VEC_DIS_NEW[r_out] = k

     OUT = GENERATE_PARA_LIST_FC(VEC_DIS_NEW,data_con,T_0,v,mu)
     MEANS_LIST = OUT[[1]]
     T_m_LIST   = OUT[[2]]
     m_LIST     = OUT[[3]]
     out = COMPUTE_BGE_FC(T_m_LIST, m_LIST, incidence, a, v, T_0)
     log_ML     = out[[1]] 

     SCORES_VEC[k] = log(m-K)-log(K) + log_ML + log(dpois(K, 1))} 
 
     VEC_DIS_NEW[r_out] = K+1
     OUT = GENERATE_PARA_LIST_FC(VEC_DIS_NEW,data_con,T_0,v,mu)
     MEANS_LIST = OUT[[1]]
     T_m_LIST   = OUT[[2]]
     m_LIST     = OUT[[3]]
     out        = COMPUTE_BGE_FC(T_m_LIST, m_LIST, incidence, a, v, T_0)
     log_ML     = out[[1]]

     SCORES_VEC[K+1] = log(K) + log_ML + log(dpois(K+1, 1))

     cons = max(SCORES_VEC) 

     SCORES_VEC = exp(SCORES_VEC - cons)

     x = sample(1:(K+1),1,replace=TRUE,prob=SCORES_VEC/sum(SCORES_VEC))
 

    VEC_DIS_NEW[r_out] = x     
    VEC_DIS            = VEC_DIS_NEW

    OUT = GENERATE_PARA_LIST_FC(VEC_DIS,data_con,T_0,v,mu)

    MEANS_LIST = OUT[[1]]
    T_m_LIST   = OUT[[2]]
    m_LIST     = OUT[[3]]
  
    ################################################################################

    out = COMPUTE_BGE_FC(T_m_LIST, m_LIST, incidence, a, v, T_0)
 
     bge_old     = out[[1]]
     P_local_num = out[[2]]
     P_local_den = out[[3]]

    ########################################################################## 

    L1[[z]] <- incidence
    L2[[z]] <- VEC_DIS
    L3[[z]] <- bge_old + compute_log_prior_allocation(VEC_DIS) + log(dpois(max(VEC_DIS), 1))

 
    }
         
  }
  return(list(L1,L2,L3))
}



##################################################################################################


compute_log_prior_allocation <- function(VEC_DIS) {


K = max(VEC_DIS)
N = length(VEC_DIS)

log_prior = 0

if (K>1) { log_prior = log(gamma(K))-log(gamma(N+K))  + sum(log(gamma(table(VEC_DIS) + 1)))  }

return(log_prior)

}



