
################################################################################

child <- function(edges,n){         # input: the numbers of the edges in the incidence matrix and the number of nodes
  p <- ceiling(edges/n)
  return(p)
}

parent <- function(edges,n){
  ch <- edges + n - child(edges,n)*n
  return(ch)
}

################################################################################

top_order <- function(incidence){
  n <- nrow(incidence)
  Order <- numeric(n)
  fan_in <- numeric(n)
  no_fan_in <- numeric(0)
  m <- 1
  for (p in 1:n){                                       # number of parent nodes at the beginning
    fan_in[p] <- sum(incidence[,p])
  }
  no_fan_in <- which(fan_in==0)
  while (length(which(Order==0))>0){                    # as long as there is a node without an order
    fan_in[which(incidence[no_fan_in[1],]==1)] <- fan_in[which(incidence[no_fan_in[1],]==1)] - 1
    no_fan_in <- c(no_fan_in, c(which(incidence[no_fan_in[1],]==1),which(fan_in==0))[duplicated(c(which(incidence[no_fan_in[1],]==1),which(fan_in==0)))])
    Order[m] <- no_fan_in[1]
    no_fan_in <- no_fan_in[-1]
    m <- m+1
  }
  return(Order)
}

################################################################################

order.edges <- function(incidence){
  top.order <- top_order(incidence)
  n <- length(top.order)
  edges <- which(incidence!=0)
  children <- child(edges,n)
  parents <- parent(edges,n)
  m <- length(edges)
  ordered_edges  <- numeric(m)
  incidence_n <- incidence
  tog <- matrix(c(edges,parents,children,ordered_edges),ncol=4, byrow=FALSE)
  k <- 1
  while(any(tog[,4]==0)){
    node1 <- top.order[which(colSums(incidence_n[,top.order])>0)][1]    # first node in top. order that has at least one parent
    par1<- tog[which(tog[,3]==node1),2]                # find the parents of  first child in the top. order that has an unordered edge incident into it
    g <- par1[which(par1>0)]
    f1 <- numeric(length(g))
    for (i in 1:length(g)){
      f1[i] <- which(top.order==g[i])
    }
    par2 <- g[which.max(f1)]                           # find the highest ordered node that has an edge leading into node1
    tog[which(tog[,2]==par2 & tog[,3]==node1),4] <- k
    k <- k + 1
    incidence_n[tog[which(tog[,2]==par2 & tog[,3]==node1),1]] <- 0     # delete the edge in the "incidence" matrix
    tog[which(tog[,2]==par2 & tog[,3]==node1),2] <- 0
  }
  to <- matrix(c(edges,parents,children,tog[,4]),ncol=4,byrow=FALSE)
  return(to)                                          # return the whole matrix, the order is the fourth column
}

#################################################################################
### DAG-to-CPDAG algorithm
#   +1 if the edge is "compelled"
#   -1 if the edge is "reversible"
#######################################################################
cpdag <- function(incidence){
  z <- order.edges(incidence)
  new_mat <- cbind(z,numeric(nrow(z)))    # edges, parents, children, order, zeros
  n_mat <- new_mat[order(new_mat[,4]),]   # sort the edges by its order
  vec <- numeric(nrow(z))
  while(any(vec==0)){                                  # while there are unlabeled edges            l.3
    if (length(vec)>1){                                  # if there are at least 2 edges
      first <- which(n_mat[,5]==0)[1]                    # first EDGE that ist labeled "unknown" (0)  l.4
      parent1 <- n_mat[first,2]                          # x   parent NODE
      child1 <- n_mat[first,3]                           # y   child NODE
      comp1 <- n_mat[which(n_mat[,3]==parent1 & n_mat[,5]==1),2]      # w NODES that have an edge incident into the parent labeled compelled)
    }
    if (length(vec)==1){
      first <- which(n_mat[5]==0)                      # first edge that ist labeled "unknown" (0)
      parent1 <- n_mat[2]                             # x   parent
      child1 <- n_mat[3]                              # y   child
      comp1 <- numeric(0)
    }
    for (j in comp1){                                   #                                            l.5
      if (incidence[j,child1]==0){                     # if w is not a parent of the child          l.6
        n_mat[first,5] <- 1                             # label x -> y compelled                     l.7
        n_mat[which(n_mat[,3]==child1),5] <- 1          # label every edge incident into y compelled l.7
        vec[first] <- 1
        vec[which(n_mat[,3]==child1)] <- 1
        break
      }
      if (incidence[j,child1]!=0)    {
        n_mat[which(n_mat[,2]==j & n_mat[,3]==child1),5] <- 1  # label w -> y compelled                l.10
        vec[which(n_mat[,2]==j & n_mat[,3]==child1)] <- 1
      }
    }
    if (length(vec)>1){
      if(n_mat[first,5]==0){
        
        moep <- n_mat[which(n_mat[,3]==child1 & n_mat[,2]!=parent1),2]      # other parents of the child
        if(length(moep)>0){                              #                     l.11
          for(o in moep){
            if(incidence[o,parent1]==0){
              vec[first] <- 1
              vec[which(n_mat[,3]==child1 & n_mat[,5]==0)] <- 1
              n_mat[first,5] <- 1                                     # label x -> y compelled
              n_mat[which(n_mat[,3]==child1 & n_mat[,5]==0),5] <- 1   # label all "unknown" edges incident into y compelled
              break
            }
            if(all(incidence[moep,parent1]!=0)){
              vec[first] <- -1
              vec[which(n_mat[,3]==child1 & n_mat[,5]==0)] <- -1
              n_mat[first,5] <- -1                                    # label x -> y reversible
              n_mat[which(n_mat[,3]==child1 & n_mat[,5]==0),5] <- -1  # label all "unknown" edges incident into y reversible
            }
          }
        }
        if(length(moep)==0){
          vec[first] <- -1
          vec[which(n_mat[,3]==child1 & n_mat[,5]==0)] <- -1
          n_mat[first,5] <- -1                                    # label x -> y reversible
          n_mat[which(n_mat[,3]==child1 & n_mat[,5]==0),5] <- -1  # label all "unknown" edges incident into y reversible
        }
      }
    }
    if (length(vec)==1){
      n_mat[5] <- -1                                    # label x -> y reversible
      vec <- -1
    }
  }
  return(n_mat)
}

##################################################################################################

EXTRACT_PRECISION_OF_DAG_NEW <- function(SIGMA,DAG){
  
  n = dim(SIGMA)[1] 
  
  V_VEC = matrix(0,n,1)
  B_MAT = matrix(0,n,n)
  
  for (i in 1:n)
  {
    parents_i = which(DAG[,i]==1)
    
    if (length(parents_i)>0)
    {
      SIGMA_parents     = SIGMA[parents_i,parents_i]
      
      
      R                 = chol(SIGMA_parents)
      R_inv             = solve(R,diag(length(parents_i)))
      
      SIGMA_parents_inv = R_inv %*% t(R_inv)
      
      
      B_MAT[parents_i,i] =             t(SIGMA[i,parents_i] %*% SIGMA_parents_inv)
      V_VEC[i,1]         = SIGMA[i,i] - SIGMA[i,parents_i] %*% SIGMA_parents_inv %*% SIGMA[parents_i,i]
    }
    else
    {V_VEC[i,1] = SIGMA[i,i]}
    
  }
  
  order_DAG = top_order(DAG)
  
  indicis = sort(order_DAG,index.return=TRUE) # To be able to re-order to the original variable order 
  
  indicis = indicis$ix
  
  V_VEC = V_VEC[order_DAG]
  B_MAT = B_MAT[order_DAG,order_DAG]
  
  W = COMPUTE_PRECISION_MATRIX_NEW(V_VEC,B_MAT)
  
  W         = W[indicis,indicis]
  R         = chol(W)
  R_inv     = solve(R,diag(n))   
  SIGMA_NEW = R_inv %*% t(R_inv)
  
  return(SIGMA_NEW) 
}

#########################################################

COMPUTE_PRECISION_MATRIX_NEW <- function(V_VEC,B_MAT)
  
{
  n = length(V_VEC)
  
  W = 1/V_VEC[1]
  
  for (i in 1:(n-1)){
    
    W = W + (B_MAT[1:i,i+1] %*% t(B_MAT[1:i,i+1]))/V_VEC[i+1]
    
    new_vec = (-1) * B_MAT[1:i,i+1]/V_VEC[i+1]
    
    W = rbind(cbind(W,new_vec),cbind(t(new_vec),1/V_VEC[i+1]))
  }
  
  return(W)
  
}

#######################################################################


top_order <- function(incidence){
  
  n <- nrow(incidence)
  Order <- numeric(n)
  fan_in <- numeric(n)
  no_fan_in <- numeric(0)
  m <- 1
  for (p in 1:n){                                       # number of parent nodes at the beginning
    fan_in[p] <- sum(incidence[,p])
  }
  no_fan_in <- which(fan_in==0)
  while (length(which(Order==0))>0){                    # as long as there is a node without an order
    fan_in[which(incidence[no_fan_in[1],]==1)] <- fan_in[which(incidence[no_fan_in[1],]==1)] - 1
    no_fan_in <- c(no_fan_in, c(which(incidence[no_fan_in[1],]==1),which(fan_in==0))[duplicated(c(which(incidence[no_fan_in[1],]==1),which(fan_in==0)))])
    Order[m] <- no_fan_in[1]
    no_fan_in <- no_fan_in[-1]
    m <- m+1
  }
  return(Order)
}

### calculation of the first ancestor matrix:
ancestor <- function(incidence){
  incidence1 <- incidence
  incidence2 <- incidence
  k <- 1
  while (k < nrow(incidence)){
    incidence1 <- incidence1%*%incidence
    incidence2 <- incidence2 + incidence1
    k <-k+1
  }
  incidence2[which(incidence2[,]>0)] <- 1
  return(t(incidence2))}


################################################################################
### function for the computation of c(n, alpha)
c_function <- function(N,A){
  fact <- numeric(N)
  for (i in 1:N){
    fact[i] <- -lgamma((A+1-i)/2)
  }
  product <- sum(fact) -(A*N/2)*log(2)- (N*(N-1)/4)*log(pi)
  return(product)}



################################################################################
### assign the topological order of the descendants of the child
des_top_order <- function(incidence, ancest1,child){
  n <- nrow(incidence)
  top <- top_order(incidence)
  position_child <- which(top==child)
  top_all_after <- top[position_child:n]                # top. order without the "first" nodes
  desc <- which(ancest1[,child]==1)                     # descendants of the child
  inter_step <- c(child,desc,top_all_after)
  des_top <- inter_step[which(duplicated(inter_step))]
  return(des_top)
}




################################################################################

cpdag_list <- function(list.inc,E){    # E: end of burnIn phase
  L <- list()
  G <- list()
  nodes <- dim(list.inc[[1]])[1]
  mat.sum <- matrix(numeric(nodes*nodes),nrow=nodes)
  for (i in E:length(list.inc)){
    k <- cpdag(list.inc[[i]])
    dummy <- matrix(numeric(nodes*nodes),nrow=nodes)
    if(length(nrow(k))!=0){
      dummy[k[,1]] <- k[,5]
      L[[i]] <- dummy
    }
    if(length(nrow(k))==0 && length(k)>0){
      dummy[k[1]] <- k[5]
      L[[i]] <- dummy
    }
    mat.com <-matrix(numeric(nodes*nodes),nrow=nodes)
    mat.re <- matrix(numeric(nodes*nodes),nrow=nodes)
    com <- which(L[[i]]>0)
    re <- which(L[[i]]<0)
    mat.com[com] <- 1
    mat.re[re] <- 1
    mat <- mat.com + mat.re + t(mat.re)
    G[[i]] <- mat
    mat.sum <- mat.sum + mat
  }
  return(list(L,G, (mat.sum/(length(list.inc)- E+1))))
}

###############################################################################################

extract_cpdag_from_dag <- function(true_incidence){    
  
  L <- list()
  
  nodes <- dim(true_incidence)[1]
  
  k <- cpdag(true_incidence)
  
  
  dummy <- matrix(numeric(nodes*nodes),nrow=nodes)
  
  if(length(nrow(k))!=0){
    dummy[k[,1]] <- k[,5]
    L <- dummy
  }
  
  if(length(nrow(k))==0 && length(k)>0){
    dummy[k[1]] <- k[5]
    L <- dummy
  }
  
  mat.com <- matrix(numeric(nodes*nodes),nrow=nodes)
  mat.re  <- matrix(numeric(nodes*nodes),nrow=nodes)
  
  com <- which(L>0)
  re  <- which(L<0)
  
  mat.com[com] <- 1
  mat.re[re] <- 1
  mat <- mat.com + mat.re + t(mat.re)
  
  return(mat)
}

#########################################################################################################





PROPOSE_NEW_DAG <- function(incidence,ancestMat,fan.in){

n = nrow(incidence)
  
### 1.) number of neighbour graphs by edge deletions
num_deletion <- sum(incidence)

### 2.) number of neighbour graphs  by edge additions    1- E(i,j) - I(i,j) - A(i,j)
inter_add <- which(matrix(rep(1,n*n),nrow=n) - diag(1,n,n) - incidence - ancestMat >0)
add <- matrix(numeric(n*n),nrow=n)
add[inter_add] <- 1
add[,which(colSums(incidence)>fan.in-1)] <- 0
num_addition <- sum(add)

### 3.) number of neighbour graphs by edge reversals    I - (I^t * A)^t
inter_rev <- which(incidence - t(t(incidence)%*% ancestMat)==1)
re <- matrix(numeric(n*n),nrow=n)
re[inter_rev] <- 1
re[which(colSums(incidence)>fan.in-1),] <- 0 
num_reversal <- sum(re)

##### total number of neighbour graphs:
total <- sum(num_deletion,num_addition,num_reversal)

### proposal probability:
p_forward <- 1/total

############## sampling a new graph (or rather sampling an edge to shift)
### sample one of the three single edge operations
random <- sample(1:total,1)

operation <- 0 # if the single edge operation is (will be) an edge reversal
if (random > total - num_reversal){operation <- 1}

#### new incidence matrix
incidence_new <- incidence

if (random <= num_deletion){             # if edge deletion 
  if(length(which(incidence>0))>1){
    new_edge <- sample(which(incidence>0),1)} # sample one edge
  else{new_edge <- which(incidence>0)}
  incidence_new[new_edge] <- 0}            # and delete it

if (random > (total - num_reversal)){      # if edge reversal was sampled
  if(num_reversal>1){
    new_edge <- sample(which(re==1),1)}  # sample one of the existing edges where a reversal leads to a valid graph
  else{
    new_edge <- which(re==1)}
  incidence_new[new_edge] <- 0             # delete it
  junk <- matrix(numeric(n*n),nrow=n)      # creating a matrix with all entries zero
  junk[new_edge] <- 1                      # an only a "1" at the entry of the new (reversed) edge
  incidence_new <- incidence_new + t(junk)}# sum the deleted matrix and the "junk-matrix"

if (random <= (total - num_reversal) & random > num_deletion){     # if edge addition was sampled
  if(num_addition>1){
    new_edge <- sample(which(add==1),1)} # sample one of the existing edges where a addition leads to a valid graph
  else{
    new_edge <- which(add==1)}
  incidence_new[new_edge] <- 1             # and add it
}


#######################################################################################################################
#################### Updating the ancestor matrix

# creating a matrix with dimensions of the incidence matrix and all entries zero except for the entry of the chosen edge
help_matrix <- matrix(numeric(n*n),nrow=n)
help_matrix[new_edge] <- 1

# numbers of the nodes that belong to the shifted edge
parent <- which(rowSums(help_matrix)==1)
child  <- which(colSums(help_matrix)==1)

### updating the ancestor matrix (after edge reversal)
## edge deletion
ancestor_new <- ancestMat
if (operation==1){
  ancestor_new[c(child,which(ancestMat[,child]==1)),] <- 0           # delete all ancestors of the child and its descendants                                           #
  top_name <- des_top_order(incidence_new, ancestMat, child)
  for (d in top_name){
    for(g in which(incidence_new[,d]==1)) {
      ancestor_new[d,c(g,(which(ancestor_new[g,]==1)))] <- 1
    }
  }
  ## edge addition
  anc_parent <- which(ancestor_new[child,]==1)                     # ancestors of the new parent
  des_child <- which(ancestor_new[,parent]==1)                     # descendants of the child
  ancestor_new[c(parent,des_child),c(child,anc_parent)] <- 1
}

### updating the ancestor matrix (after edge deletion)
if (random <= num_deletion){
  ancestor_new[c(child,which(ancestMat[,child]==1)),] <- 0           # delete all ancestors of the child and its descendants                                           #
  top_name <- des_top_order(incidence_new, ancestMat, child)
  for (d in top_name){
    for(g in which(incidence_new[,d]==1)) {
      ancestor_new[d,c(g,(which(ancestor_new[g,]==1)))] <- 1
    }
  }
}

# updating the ancestor matrix (after edge addition)
if (random <= total - num_reversal & random > num_deletion){
  anc_parent <- which(ancestMat[parent,]==1)        # ancestors of the new parent
  des_child <- which(ancestMat[,child]==1)          # descendants of the child
  ancestor_new[c(child,des_child),c(parent,anc_parent)] <- 1
}




######################################################################################################
####### ... the number of neighbour graphs/proposal probability for the proposed graph
### 1.) number of neighbour graphs obtained by edge deletions
num_deletion_new <- sum(incidence_new)

### number of neighbour graphs obtained by edge additions    1- E(i,j) - I(i,j) - A(i,j)
inter_add.new <- which(matrix(rep(1,n*n),nrow=n) - diag(1,n,n) - incidence_new - ancestor_new >0)
add.new <- matrix(numeric(n*n),nrow=n)
add.new[inter_add.new] <- 1
add.new[,which(colSums(incidence_new)>fan.in-1)] <- 0
num_addition_new <- sum(add.new)

### number of neighbour graphs obtained by edge reversals    I - (I^t * A)^t
inter_rev.new <- which(incidence_new - t(t(incidence_new)%*% ancestor_new)==1)
re.new <- matrix(numeric(n*n),nrow=n)
re.new[inter_rev.new] <- 1
re.new[which(colSums(incidence_new)>fan.in-1),] <- 0  # CORRECTED!!!???!!!
num_reversal_new <- sum(re.new)

##### total number of neighbour graphs:
total_new <- sum(num_deletion_new,num_addition_new,num_reversal_new)

### proposal probability:
p_backward <- 1/total_new

return(list(incidence_new,ancestor_new,p_forward,p_backward,operation,child,parent))
}



################################################################################################################



co_allocation_matrix <- function(out_mGBN){


n_samples = length(out_mGBN[[2]])

VEC_DIS_1 = out_mGBN[[2]][[1]]

n_obs = length(VEC_DIS_1)


MATRIX = matrix(0,n_obs,n_obs)

Counter = 0

for (i_sample in round((n_samples/2)+1):n_samples){


VEC_DIS = out_mGBN[[2]][[i_sample]]

for (i in 1:n_obs) {
for (j in 1:n_obs){
if(VEC_DIS[i]==VEC_DIS[j]) {MATRIX[i,j] = MATRIX[i,j] + 1}
}
}

Counter = Counter + 1

}


MATRIX = MATRIX/Counter

return(MATRIX)

}
 















