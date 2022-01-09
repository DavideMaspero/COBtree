# remove any indistinguishable variant from input data
check.indistinguishable <- function( data ) {
  
  # check for indistinguishable events
  indistinguishable <- as.numeric(which(duplicated(t(data))))
  if(length(indistinguishable)>0) {
    # merge names of indistinguishable events
    valid_colnames <- colnames(data)[-indistinguishable]
    new_colnames <- valid_colnames
    invalid_colnames <- colnames(data)[indistinguishable]
    for(i in invalid_colnames) {
      for(j in valid_colnames) {
        if(all(is.na(data[,i])==is.na(data[,j]))) { # if NAs in i and j are the same
          if(all(data[,i]==data[,j],na.rm=TRUE)) { # if not-NAs entries in i and j are the same
            new_colnames[which(valid_colnames==j)] <- paste0(new_colnames[which(valid_colnames==j)],"|",invalid_colnames[which(invalid_colnames==i)])
            next;
          }
        }
      }
    }
    # remove indistinguishable events from data
    data <- data[,valid_colnames,drop=FALSE]
    colnames(data) <- new_colnames
  }
  
  # return data
  return(data)
  
}

# build a phylogenetic tree from a variants tree (adjacency_matrix) and related samples attachments (samples_attachments)
get.phylo <- function( adjacency_matrix, valid_genotypes, samples_attachments ) {
  
  # compute Manhattan distance among valid genotypes
  distance_genotypes <- as.matrix(dist(valid_genotypes,method="manhattan"))
  
  # consider the variants tree (adjacency_matrix) to build inner nodes structure of VERSO phylogenetic tree
  edges <- which(adjacency_matrix==1,arr.ind=TRUE)
  parents_list <- rownames(adjacency_matrix)[as.numeric(edges[,"row"])]
  children_list <- rownames(adjacency_matrix)[as.numeric(edges[,"col"])]
  edges <- cbind(parents_list,children_list)
  edges_weights <- NULL
  for(i in seq_len(nrow(edges))) {
    curr_p <- rownames(valid_genotypes)[which(colnames(valid_genotypes)==edges[i,1])]
    curr_c <- rownames(valid_genotypes)[which(colnames(valid_genotypes)==edges[i,2])]
    edges_weights <- c(edges_weights,distance_genotypes[curr_p,curr_c])
  }
  nodes_list <- unique(as.vector(edges))
  Nnode <- length(nodes_list)
  parents_list_numeric <- parents_list
  children_list_numeric <- children_list
  for(nlist in seq_len(Nnode)) {
    curr_nl <- which(parents_list==nodes_list[nlist])
    if(length(curr_nl)>0) {
      parents_list_numeric[curr_nl] <- nlist
    }
    curr_nl <- which(children_list==nodes_list[nlist])
    if(length(curr_nl)>0) {
      children_list_numeric[curr_nl] <- nlist
    }
  }
  edges <- cbind(as.numeric(parents_list_numeric),as.numeric(children_list_numeric))
  colnames(edges) <- c("Parent","Child")
  
  # now consider samples attachments
  attachments <- samples_attachments[,"Genotype"]
  tip.label <- c(names(attachments),"Reference_Genotype")
  overhead <- length(tip.label)
  edges <- edges + overhead
  for(i in seq_len(length(attachments))) {
    curr_att <- which(nodes_list==colnames(valid_genotypes)[which(rownames(valid_genotypes)==attachments[i])])+overhead
    edges <- rbind(edges,c(curr_att,i))
    edges_weights <- c(edges_weights,0)
  }
  edges <- rbind(edges,c(as.numeric(edges[1,"Parent"]),(i+1)))
  edges_weights <- c(edges_weights,0)
  rownames(edges) <- seq_len(nrow(edges))
  
  # build VERSO phylogenetic tree
  phylogenetic_tree <- list(edge=edges,tip.label=tip.label,Nnode=Nnode,edge.length=edges_weights,node.label=nodes_list,root.edge=0)
  class(phylogenetic_tree) <- "phylo"
  phylogenetic_tree <- ape::keep.tip(phylogenetic_tree,tip.label)
  
  # return VERSO phylogenetic tree
  return(phylogenetic_tree)
  
}

# convert B to an adjacency matrix
as.adj.matrix <- function( B, sorting = FALSE) {
  
  # create the data structure where to save the adjacency matrix obtained from B
  adj_matrix <- array(0L,dim(B))
  rownames(adj_matrix) <- colnames(B)
  colnames(adj_matrix) <- colnames(B)
  
  # set arcs in the adjacency matrix
  for(i in seq_len((nrow(B)-1))) {
    for(j in ((i+1):nrow(B))) {
      if(all(B[i,seq_len(i)]==B[j,seq_len(i)])&&(sum(B[i,])==(sum(B[j,])-1))) {
        adj_matrix[i,j] <- 1
      }
    }
  }
  
  if(sorting == TRUE) {
    #keeping root on the first place
    ord <- c(1,1+order(colnames(adj_matrix)[-1]))
    adj_matrix <- adj_matrix[ord,ord]
  }
  
  # return the adjacency matrix obtained from B
  return(adj_matrix)
  
}

# build B from an adjacency matrix where we assume genotypes and mutations to be both ordered
as.B <- function( adj_matrix, D ) {
  
  # build data structure to save results
  n_clones <- nrow(adj_matrix)
  B <- diag(n_clones)
  
  # build B
  for(k in seq_len(n_clones)) {
    idx_child <- which(adj_matrix[k,]==1,arr.ind=TRUE)
    if(length(idx_child)==1) {
      B[idx_child,] <- B[k,] + B[idx_child,]
    }
    else if(length(idx_child)>1) {
      B[idx_child,] <- sweep(B[idx_child,],2,B[k,],"+")
    }
  }
  rownames(B) <- c("r",seq_len((nrow(B)-1)))
  mycolnames <- "r"
  for(i in 2:nrow(adj_matrix)) {
    mycolnames <- c(mycolnames,as.character(which(rownames(adj_matrix)[i]==colnames(D))))
  }
  colnames(B) <- mycolnames
  
  # return B
  return(B)
  
}


# New version
as.B.rooted <- function(rooted_adjM) {
  
  #get par-child pairs
  idx_edges <- adj.matrix.2.vect(adjM = rooted_adjM)
  
  nMut = nrow(rooted_adjM)
  
  mut_par = c()
  mut_child = mut_par
  
  for(idx in idx_edges) {
    
    r_i = ((idx-1) %% nrow(rooted_adjM)) + 1
    c_i = floor((idx-1) / nrow(rooted_adjM)) + 1
    
    mut_par <- c(mut_par, rownames(rooted_adjM)[r_i])
    mut_child <- c(mut_child, colnames(rooted_adjM)[c_i])
    
  }
  
  B <- array(data=0L,dim =c(nMut,nMut)) 
  
  colnames(B) <- colnames(rooted_adjM)
  
  #colnames(unrootB) <- as.character(vTree[1:nMut])
  #edges <- as.numeric(vTree[(nMut+1):(nMut*2)])
  
  #B must have diagonal = 1
  diag(B) <- 1
  
  next_node = unique(mut_par[!(mut_par %in% mut_child)])
  
  if(length(next_node) != 1) {
    stop("Found none or multiple roots in the adjacency matrix")
  }
  
  i = 1
  
  idx_P <- which(colnames(B)==next_node)
  colnames(B)[c(1,idx_P)] <- colnames(B)[c(idx_P,1)]
  
  
  while(i < (nMut)) {
    
    tmp_nn <- c()
    
    for(nn in next_node) {
      
      idx_C = which(mut_par==nn)
      
      if(length(idx_C)==0) {next}
      
      idx_P <- which(colnames(B)==nn)
      
      nC <- length(idx_C)
      
      B[(i+1):(i+nC),1:idx_P] <- rep(B[idx_P,1:idx_P], each = nC)
      
      colnames(B)[colnames(B) %in% mut_child[idx_C]] <- NA
      
      colnames(B)[(i+1):(i+nC)] <- mut_child[idx_C]
      
      tmp_nn <- c(tmp_nn, mut_child[idx_C])
      
      i = i + nC
      
    }
    
    next_node <- tmp_nn
    
  }
  
  #colnames(B)[1] <- 'r'
  
  return(B)
  
}


# draw B
draw.B <- function( B, mut_label = colnames(B)[-1], last_mut_node_label =  TRUE) {
  
  Broot <- Node$new('r')
  Broot$mut <- B[1,]
  
  nClone <- nrow(B)
  Clones <- list(Broot)
  
  for(rP in 1:(nrow(B)-1)) {
    
    for(rC in ((rP+1):nrow(B))) {
      
      if(all(Clones[[rP]]$mut[1:rP]==B[rC,1:rP])&&(sum(Clones[[rP]]$mut)==(sum(B[rC,])-1))) {
        if(last_mut_node_label) {
          mutName <- tail(mut_label[which(B[rC,-1]==1)], 1)
        } else {
          mutName <- paste(mut_label[which(B[rC,-1]==1)],collapse="")  
        }
        
        Clones[[rC]] <- Clones[[rP]]$AddChild(mutName)
        Clones[[rC]]$mut <- B[rC,]
        
      }
      
    }
    
  }
  
  return(Broot)
  
}

adj.matrix.2.vect <- function(adjM) {
  vect <- which(adjM != 0)
  storage.mode(vect) <- "integer"
  return(vect)
}

vect.2.adj.matrix <- function(vect = NULL, value = 1L, mut_names = NULL) {
  #IMPORTANT the matrix must be ordered
  adjM <- matrix(data = 0L, ncol = (length(vect)+1), nrow = (length(vect)+1))
  adjM[vect] <- value
  
  if(!is.null(mut_names)) {
    
    colnames(adjM) <- mut_names
    rownames(adjM) <- mut_names
    
  }
  
  return(adjM)  
  
}

edge.table.2.adj.matrix <- function(edge_table = NULL, mut_names=NULL) {
  adjM <- matrix(data = 0L, ncol = (length(mut_names)), nrow = (length(mut_names)))
  
  for(i in 1:length(edge_table)) {
    edge_count_i <- edge_table[i]
    idx_i <- as.numeric(names(edge_count_i))
    adjM[idx_i] <- as.numeric(edge_count_i)
  }
  
  if(!is.null(mut_names)) {
    
    colnames(adjM) <- mut_names
    rownames(adjM) <- mut_names
    
  }
  
  return(adjM)
  
}


linMap <- function(x, from, to) {
  mapped_x <- (((x - min(x)) ) / max(x - min(x)) * (to - from)) + from
  return(mapped_x)
}

getG <- function(inference=NULL) {
  G <- matrix(data = 0, 
              nrow = length(inference$C), 
              ncol = ncol(inference$B)-1, 
              dimnames = list(rownames(inference$C), 
                              colnames(inference$B)[-1])
  )
  
  for(i in 1:length(inference$C)) {
    
    cell_i <- inference$C[i]
    
    G[i,] <- inference$B[cell_i,2:ncol(inference$B)]
    
  }
  
  G <- G[,sort(colnames(G))]
  return(G)
}

compareD_G <- function(D=D,G=G) {
  
  D <- D[order(rownames(D)), 
         order(colnames(D))]
  D[which(D < 0)] <- NA
  
  G <- G[order(rownames(G)), 
         order(colnames(G))]
  
  diffD_G <- D - G
  
  FP = sum(diffD_G == 1, na.rm = T)
  FN = sum(diffD_G == -1, na.rm = T)
  
  FPrate = FP / sum(G == 0 & !is.na(D))
  FNrate = FN / sum(G == 1 & !is.na(D))
  
  return(list(FPrate = FPrate, FNrate = FNrate))
}

scite.2.adjM <- function(scite_ml0_fn = NULL, mutNames = NULL) {
  
  #parsing scite.out.gv
  
  #gv <- readLines(paste0(outdir,"/scite.out_ml0.gv"))
  gv <- readLines(scite_ml0_fn)
  gv <- gv[3:(length(gv)-1)]
  
  node_pairs <- t(sapply(X = gv, 
                         FUN = function(x) 
                         { 
                           x <- gsub(x=x, pattern=";",replacement="")
                           pairs <- as.numeric(unlist(strsplit(x = x, split = " -> ")))
                           return(pairs)
                         }, 
                         USE.NAMES = F, simplify = T)
  )
  
  root = max(node_pairs)
  
  #replace label
  node_pairs[node_pairs==root] <- 0
  
  nMut = max(node_pairs)
  
  #generate scite_adj
  scite_adj <- matrix(0, nMut+1, nMut+1)
  scite_adj[node_pairs+1] <- 1
  
  colnames(scite_adj) <- c(0:nMut)
  
  if(!is.null(mutNames)) {
    colnames(scite_adj) <- c('r',mutNames)[as.numeric(colnames(scite_adj))+1]
    
    mut_order <- c('r',sort(mutNames))
    idx_ordered <- match(mut_order, colnames(scite_adj))
    
    scite_adj <- scite_adj[idx_ordered,idx_ordered]
    
  } else {
    scite_adj <- scite_adj[order(colnames(scite_adj)),order(rownames(scite_adj))]
    colnames(scite_adj)[1] <- rownames(scite_adj)[1] <- 'r'
  }
  
  rownames(scite_adj) <- colnames(scite_adj)
  
  return(scite_adj)
  
}

sasc.2.adjM <- function(sasc_mlt_fn = NULL) {
  
  mlt_f <- readLines(sasc_mlt_fn)
  mlt <- strsplit(x = mlt_f, split = "\t", fixed = T)[[1]][2]
  parent_child_pairs <- strsplit(x = mlt, split = ';')[[1]]
  
  sasc_adjM <- matrix(data = 0, nrow = 1, ncol = 1, dimnames = list("germline","germline"))
  
  for(pc_i in parent_child_pairs) {
    parent_child_id <- strsplit(x = pc_i, split = '>', fixed = T)[[1]]
    
    # check if the node is already included in the adjacency matrix otherwise add an empty columns
    if(!(parent_child_id[1] %in% colnames(sasc_adjM))) {
      sasc_adjM <- cbind(sasc_adjM, rep(0, nrow(sasc_adjM)))
      sasc_adjM <- rbind(sasc_adjM, rep(0, ncol(sasc_adjM)))
      
      colnames(sasc_adjM)[ncol(sasc_adjM)] <- rownames(sasc_adjM)[nrow(sasc_adjM)] <- parent_child_id[1]
    }
    
    if(!(parent_child_id[2] %in% colnames(sasc_adjM))) {
      sasc_adjM <- cbind(sasc_adjM, rep(0, nrow(sasc_adjM)))
      sasc_adjM <- rbind(sasc_adjM, rep(0, ncol(sasc_adjM)))
      
      colnames(sasc_adjM)[ncol(sasc_adjM)] <- rownames(sasc_adjM)[nrow(sasc_adjM)] <- parent_child_id[2]
    }
    
    sasc_adjM[parent_child_id[1],parent_child_id[2]] <- 1
    
  }
  
  colnames(sasc_adjM)[which(colnames(sasc_adjM)=="germline")] <- 'r'
  rownames(sasc_adjM)[which(rownames(sasc_adjM)=="germline")] <- 'r'
  
  return(sasc_adjM)
  
}

get.ancestor.nodes <- function(adjM = NULL, t_node = NULL) {
  
  ancestors <- c()
  
  n_node = rownames(adjM)[which(adjM[,t_node] == 1)]
  
  while(length(n_node) > 0 && n_node != 'r') {
    ancestors <- c(ancestors, n_node)
    n_node <- rownames(adjM)[which(adjM[,n_node] == 1)]
  }
  return(ancestors)
}

extend.adj.matrices<-function(adjA = NULL, adjB = NULL){
  # columns not included in adjB
  mut_not_present <- setdiff(colnames(adjA),colnames(adjB))
  for(mut in mut_not_present){
    adjB <- cbind(adjB, rep(0, nrow(adjB)))
    adjB <- rbind(adjB, rep(0, ncol(adjB)))
    colnames(adjB)[ncol(adjB)] <- rownames(adjB)[nrow(adjB)] <- mut
  }
  
  # columns not included in adjA
  mut_not_present <- setdiff(colnames(adjB),colnames(adjA))
  for(mut in mut_not_present){
    adjA <- cbind(adjA, rep(0, nrow(adjA)))
    adjA <- rbind(adjA, rep(0, ncol(adjA)))
    colnames(adjA)[ncol(adjA)] <- rownames(adjA)[nrow(adjA)] <- mut
  }
  
  #sorting
  adjA <- adjA[order(colnames(adjA)),order(rownames(adjA))]
  adjB <- adjB[order(colnames(adjB)),order(rownames(adjB))]
  
  return(list(adjA = adjA, adjB = adjB))
  
}

distance.parent_child.SCITEsamples <- function(pvA = NULL, pvB = NULL) {
  parent_child_dist = sum(pvA != pvB)
  return(parent_child_dist)
}

getClonalGenotypes_sorted <- function(B = NULL) {
  clonal_geno <- B[-1,-1]
  rownames(clonal_geno) <- colnames(clonal_geno)
  
  clonal_geno <- clonal_geno[order(rownames(clonal_geno)),
                             order(colnames(clonal_geno))]
  return(clonal_geno)
}








