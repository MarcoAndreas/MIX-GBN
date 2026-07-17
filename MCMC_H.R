MCMC_H <- function(data,incidence,iterations,step_save){
  
  data = t(data)
  

  n <- nrow(data)    # number of nodes
  m <- ncol(data)    # number of observations
  
  fan.in = n-1
  
  v=1
  
  mu=numeric(n)
  
  a=n+2 

  T_0 = diag(1,n,n)
  
  T_0 <- (a-n-1) * T_0 # scaling T_0
  
  L1 <- list()    # DAGs
  L2 <- list()    # log scores
  L3 <- list()
  L4 <- list() 

  ################################################################################
  ### computation of the log BGe Score of the FIRST graph
 
  T_m = T_0 + (m-1)* cov(t(data)) + ((v*m)/(v+m))* (mu - rowMeans(data))%*%t(mu - rowMeans(data))
  
  P_local_num = numeric(n) 
  P_local_den = numeric(n) 
    
 
  out = COMPUTE_BGE_H(T_m,m,P_local_num,P_local_den,incidence,a,v,T_0)
 
  bge_old     = out[[1]]
  P_local_num = out[[2]]
  P_local_den = out[[3]]

  L1[[1]] <- incidence                    
  L2[[1]] <- bge_old                     

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
      
      out = COMPUTE_BGE_LOCAL_H(operation,child_node,parent_node,T_m,m,P_local_num,P_local_den,incidence_new,a,v,T_0)
      
      bge_new         = out[[1]]
      P_local_num_new = out[[2]]
      P_local_den_new = out[[3]]
      
      acceptance <- exp((bge_new +log(p_backward)) - (bge_old +log(p_forward)))
      
      rand <- runif(1)
      
      if(acceptance > rand){
        incidence    <- incidence_new
        ancestMat    <- ancestMat_new
        bge_old      <- bge_new
        P_local_num  <- P_local_num_new
        P_local_den  <- P_local_den_new
      }
    }
    
    out_para = SAMPLE_PARAMETERS_H(T_m,m,data,a,v,mu,incidence)
  
    L1[[z]] <- incidence
    L2[[z]] <- bge_old
    L3[[z]] <- out_para[[1]]
    L4[[z]] <- out_para[[2]]
  }
  return(list(L1,L2,L3,L4))
  
  
}
