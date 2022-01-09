computeMSTweighted <- function(trees=NULL, likelihood=NULL, ncores = 1) {
  
  # COMPUTE CONSENSUS TREE WITH GRAPH-BASED APPROACH
  all_trees_df <- do.call(rbind,  parallel::mclapply(trees, function(x) {
    adjM <- as.adj.matrix(x, sorting = TRUE)
    sparse <- adj.matrix.2.vect(adjM)
    return(sparse)
  }, mc.cores = ncores))
  
  mut_names <- colnames(as.adj.matrix(B = trees[[1]], sorting = T))
  
  idx_uniq_sol <- which(!duplicated(all_trees_df))
  
  #creare consensus adj
  
  unique_trees <- trees[idx_uniq_sol]
  unique_lik <- likelihood[idx_uniq_sol]
  unique_trees_df <- all_trees_df[idx_uniq_sol,]
  
  #calcolare MST
  
  cons_adjM <- matrix(data = 0, 
                      nrow = length(mut_names), 
                      ncol = length(mut_names), 
                      dimnames = list(mut_names,mut_names)
  )
  
  lik_0_100 <- linMap(x = unique_lik, 
                      from = 1, 
                      to = 100)
  
  min_log_lik <- min(lik_0_100)
  for(i in seq(1, length(unique_trees))) {
    adjM <- vect.2.adj.matrix(vect = unique_trees_df[i,], 
                              value = exp(lik_0_100[i]),
                              mut_names = mut_names)
    #adjM <- as.adj.matrix(B = unique_trees[[i]], 
    #                      sorting = TRUE)
    cons_adjM <- cons_adjM + adjM
  }
  
  # Generate graph from consensus adjacency matrix
  
  g <- new("graphNEL", nodes=colnames(cons_adjM), edgemode="directed")
  
  idx_edges <- adj.matrix.2.vect(cons_adjM)
  
  for(idx in idx_edges) {
    
    r_i = ((idx-1) %% nrow(cons_adjM)) + 1
    c_i = floor((idx-1) / nrow(cons_adjM)) + 1
    
    g <- graph::addEdge(from = rownames(cons_adjM)[r_i], 
                        to = colnames(cons_adjM)[c_i], 
                        graph = g, 
                        weights = cons_adjM[idx] / max(cons_adjM))
    
  }
  
  min_span_tree <- RBGL::edmondsOptimumBranching(g)
  
  min_span_adjM  <- matrix(data = 0, 
                           nrow = length(mut_names), 
                           ncol = length(mut_names), 
                           dimnames = list(mut_names,mut_names)
  )
  
  for(i in seq(1, ncol(min_span_tree$edgeList))) {
    r_i <- as.character(min_span_tree$edgeList["from", i])
    c_i <- as.character(min_span_tree$edgeList["to", i])
    min_span_adjM[r_i, c_i] <- min_span_tree$weights[i]
    
  }
  
  cons_tree <- list()
  
  cons_tree$min_span_tree <- as.B.rooted(rooted_adjM = min_span_adjM)
  cons_tree$min_span_weigth_adj = min_span_adjM
  cons_tree$cons_adj <- cons_adjM
  
  return(cons_tree)
  
}

computeMST.fromTreeList <- function(trees=NULL, ncores = 1) {
  
  mut_names <- colnames(as.adj.matrix(B = trees[[1]], sorting = T))
  
  cons_adjM <- matrix(data = 0, 
                      nrow = length(mut_names), 
                      ncol = length(mut_names), 
                      dimnames = list(mut_names,mut_names)
  )
  
  for(tree in trees) {
    cons_adjM <- cons_adjM + as.adj.matrix(tree, sorting = TRUE)
  }
  
  
  # cons_adjM <- Reduce('+', 
  #                     parallel::mclapply(trees, function(x){as.adj.matrix(B=x, sorting=TRUE)}, mc.cores = ncores)
  #                     )
  # 
  
  #calcolare MST
  
  g <- new("graphNEL", nodes=colnames(cons_adjM), edgemode="directed")
  
  idx_edges <- adj.matrix.2.vect(cons_adjM)
  
  for(idx in idx_edges) {
    
    r_i = ((idx-1) %% nrow(cons_adjM)) + 1
    c_i = floor((idx-1) / nrow(cons_adjM)) + 1
    
    g <- graph::addEdge(from = rownames(cons_adjM)[r_i], 
                        to = colnames(cons_adjM)[c_i], 
                        graph = g, 
                        weights = cons_adjM[idx] / max(cons_adjM))
    
  }
  
  min_span_tree <- RBGL::edmondsOptimumBranching(g)
  
  min_span_adjM  <- matrix(data = 0, 
                           nrow = length(mut_names), 
                           ncol = length(mut_names), 
                           dimnames = list(mut_names,mut_names)
  )
  
  for(i in seq(1, ncol(min_span_tree$edgeList))) {
    r_i <- as.character(min_span_tree$edgeList["from", i])
    c_i <- as.character(min_span_tree$edgeList["to", i])
    min_span_adjM[r_i, c_i] <- min_span_tree$weights[i]
    
  }
  
  cons_tree <- list()
  
  cons_tree$min_span_tree <- as.B.rooted(rooted_adjM = min_span_adjM)
  cons_tree$min_span_weigth_adj = min_span_adjM
  cons_tree$cons_adj <- cons_adjM
  
  return(cons_tree)
  
}

computeMST.fromAdjDF <- function(adjDF=NULL, mut_names = NULL) {
  
  
  edge_count <- table(adjDF)
  cons_adjM <- edge.table.2.adj.matrix(edge_table = edge_count, mut_names = mut_names)
  
  #  for(tree in trees) {
  #    cons_adjM <- cons_adjM + as.adj.matrix(tree, sorting = TRUE)
  #  }
  
  #calcolare MST
  
  g <- new("graphNEL", nodes=mut_names, edgemode="directed")
  
  idx_edges <- adj.matrix.2.vect(cons_adjM)
  
  for(idx in idx_edges) {
    
    r_i = ((idx-1) %% nrow(cons_adjM)) + 1
    c_i = floor((idx-1) / nrow(cons_adjM)) + 1
    
    g <- graph::addEdge(from = rownames(cons_adjM)[r_i], 
                        to = colnames(cons_adjM)[c_i], 
                        graph = g, 
                        weights = cons_adjM[idx] / max(cons_adjM))
    
  }
  
  min_span_tree <- RBGL::edmondsOptimumBranching(g)
  
  min_span_adjM  <- matrix(data = 0, 
                           nrow = length(mut_names), 
                           ncol = length(mut_names), 
                           dimnames = list(mut_names,mut_names)
  )
  
  for(i in seq(1, ncol(min_span_tree$edgeList))) {
    r_i <- as.character(min_span_tree$edgeList["from", i])
    c_i <- as.character(min_span_tree$edgeList["to", i])
    min_span_adjM[r_i, c_i] <- min_span_tree$weights[i]
    
  }
  
  cons_tree <- list()
  
  cons_tree$min_span_tree <- as.B.rooted(rooted_adjM = min_span_adjM)
  cons_tree$min_span_weigth_adj = min_span_adjM
  cons_tree$cons_adj <- cons_adjM
  
  return(cons_tree)
  
}

compute.consensus <- function(cons_adjM = cons_adjM, tot_tree = tot_selected_trees) {
  
  g <- new("graphNEL", nodes=colnames(cons_adjM), edgemode="directed")
  
  idx_edges <- adj.matrix.2.vect(cons_adjM)
  
  for(idx in idx_edges) {
    
    r_i = ((idx-1) %% nrow(cons_adjM)) + 1
    c_i = floor((idx-1) / nrow(cons_adjM)) + 1
    
    g <- graph::addEdge(from = rownames(cons_adjM)[r_i], 
                        to = colnames(cons_adjM)[c_i], 
                        graph = g, 
                        weights = cons_adjM[idx] / tot_tree)
    
  }
  
  max_spanning_tree <- RBGL::edmondsOptimumBranching(g)
  
  max_span_adjM  <- cons_adjM*0
  
  for(i in seq(1, ncol(max_spanning_tree$edgeList))) {
    r_i <- as.character(max_spanning_tree$edgeList["from", i])
    c_i <- as.character(max_spanning_tree$edgeList["to", i])
    max_span_adjM[r_i, c_i] <- max_spanning_tree$weights[i]
    
  }
  
  cons_tree <- list()
  
  cons_tree$cob_tree <- as.B.rooted(rooted_adjM = max_span_adjM)
  cons_tree$cob_weigth_adj = max_span_adjM
  cons_tree$cob_adj <- cons_adjM
  cons_tree$sampled_trees <- tot_tree
  
  return(cons_tree)
  
}


computeMST.fromSCITEsamples <- function(scite_samples_fn = NULL, mut_names = NULL, tree_every = 1, logLik_thr = 0.2, tree_buffer = 100000, verbose=TRUE) {
  
  
  # Read scite.samples file
  if(!file.exists(scite_samples_fn)) {
    stop("scite.samples not found")
  } else {
    tree_list_raw <- readLines(scite_samples_fn)
  }
  
  logLikelihoodV <-  sapply(tree_list_raw, function(x){as.numeric(strsplit(x,'\t')[[1]][1])},USE.NAMES=F)
  minLogLikelihood <- max(logLikelihoodV)*(1+logLik_thr)
  
  # Filter out tree with log likelihood lower than minLogLikelihood
  
  idx_keep <- logLikelihoodV >= minLogLikelihood
  tree_list_raw <- tree_list_raw[idx_keep]
  logLikelihoodV <- logLikelihoodV[idx_keep]
  
  
  # Process each tree to get the consensus matrix
  
  parent_vect = strsplit(x = tree_list_raw[1], split = '\t')[[1]][3]
  nMut = length(strsplit(x = parent_vect, split = ' ')[[1]])
  
  if(is.null(mut_names)) {
    mut_names <- 1:nMut
  }
  
  if(tree_buffer >= length(tree_list_raw)) {
    i_seq <- c(0, length(tree_list_raw))
  } else {
    i_seq <- union(seq(0, length(tree_list_raw), by = tree_buffer), length(tree_list_raw)) 
  }
  
  cons_adjM <- matrix(data = 0, nrow = nMut + 1, ncol = nMut+1)
  colnames(cons_adjM) <- rownames(cons_adjM) <- c(mut_names,'r')
  
  #pb <- progress::progress_bar$new(
  #      format = "  Processing trees [:bar] :current/:total (:percent) eta ~ :eta",
  #      total = max(i_seq), clear = FALSE, width= 80)
  
  
  total_trees <- max(i_seq)
  
  if(verbose) {
    count = 0
    cat('\n', 
        'Processing trees.. ', count, 
        '/', total_trees, 
        ' (',format(count/total_trees*100, digits = 2, nsmall=2),'%)',
        sep=''
    )
  }
  for(i in 1:(length(i_seq)-1)) {
    f_tree = i_seq[i]+1
    e_tree = i_seq[i+1]
    
    tmp_parent_vector_matrix <- matrix(data = 0, nrow = (1 + e_tree - f_tree), ncol = nMut)
    
    tmp_r <- 1
    
    for(tree_i in seq(f_tree,e_tree,by=tree_every)) {
      
      tree_raw_i <- strsplit(x = tree_list_raw[tree_i], split = '\t')[[1]]
      tmp_parent_vector_matrix[tmp_r,] <- as.numeric(strsplit(x = tree_raw_i[3], split = ' ')[[1]]) + 1
      
      tmp_r = tmp_r + 1
      
      if(verbose){
        count = count + 1
        if((count %% 1000) == 0 || count == total_trees) {
          cat('\r',
              'Processing trees.. ', count, 
              '/', total_trees, 
              ' (',format(count/total_trees*100, digits = 2, nsmall = 2),'%)',
              sep=''
          )
        }
      }
      
    }
    if(verbose){
      cat('\n', 
          'Processing nodes.. ', 1, 
          '/', nMut, 
          ' (',format(1/nMut*100, digits = 2, nsmall=2),'%)',
          '          ',
          sep=''
      )
    }
    for(mut_i in 1:nMut) {
      
      parent_i_freq <- table(tmp_parent_vector_matrix[,mut_i])
      
      for(pr_i in 1:length(parent_i_freq)) {
        
        p_index <- as.numeric(names(parent_i_freq)[pr_i])
        edge_count <- as.numeric(parent_i_freq[pr_i])
        
        cons_adjM[p_index,mut_i] <- cons_adjM[p_index,mut_i] + edge_count
        
      } 
      if(verbose) {
        cat('\r', 
            'Processing nodes.. ', mut_i, 
            '/', nMut, 
            ' (',format(mut_i/nMut*100, digits = 2, nsmall=2),'%)',
            '          ',
            sep=''
        )
      }
    }
  }
  
  if(verbose) {
    cat('\n')
  }
  
  cons_tree <- compute.consensus(cons_adjM = cons_adjM,
                                 tot_tree = length(tree_list_raw))
  
  
  return(cons_tree)
  
}  

computeMST.fromSASCsamples <- function(sasc_samples_fn = NULL, tree_every = 1, logLik_thr = 0.2, tree_buffer = 10000, keep_only_ml_events = T, fix_inconsistencies = F, p_rem_bm = 0.05) {
  
  
  # Read scite.samples file
  if(!file.exists(sasc_samples_fn)) {
    stop("sasc.samples not found")
  } else {
    tree_list_raw <- readLines(sasc_samples_fn)
  }
  
  logLikelihoodV <-  unlist(lapply(tree_list_raw, function(x){strsplit(x,'\t')[[1]]}))
  
  rm(tree_list_raw)
  
  sampledTreesV <- logLikelihoodV[seq(2,length(logLikelihoodV), by=2)]
  logLikelihoodV <- as.numeric(logLikelihoodV[seq(1,length(logLikelihoodV), by=2)])
  
  minLogLikelihood_thr <- max(logLikelihoodV)*(1+logLik_thr)
  
  # Filter out tree with log likelihood lower than minLogLikelihood_thr
  
  idx_keep <- logLikelihoodV >= minLogLikelihood_thr
  sampledTreesV <- sampledTreesV[idx_keep]
  logLikelihoodV <- logLikelihoodV[idx_keep]
  
  tot_selected_trees <- length(sampledTreesV)
  
  # Count parent_child_node_pairs to build the consensus matrix
  
  if(keep_only_ml_events) {
    
    # get information about the mutation included in the ML trees
    
    ML_tree_pairsV <- unlist(lapply(sampledTreesV[which(logLikelihoodV == max(logLikelihoodV))], function(x){strsplit(x,';')[[1]]}))
    ML_all_mut <- unique(unlist(lapply(ML_tree_pairsV, 
                                       function(x){
                                         p <- strsplit(x,'>')[[1]]
                                         return(p)
                                       }
    )
    )
    )
    
    ML_back_mutations <- ML_all_mut[startsWith(x = ML_all_mut, prefix = '-1')]
    possible_bm <- setdiff(x = paste0("-", ML_all_mut[!startsWith(x = ML_all_mut, prefix = '-')]), y = c(ML_back_mutations, "-germline"))
  }
  
  cons_adjM <- matrix(data = 0, nrow = 1, ncol = 1, dimnames = list("germline","germline"))
  
  if(tree_buffer >= tot_selected_trees) {
    i_seq <- c(0, tot_selected_trees)
  } else {
    i_seq <- union(seq(0, tot_selected_trees, by = tree_buffer), tot_selected_trees) 
  }
  
  actual_sampled_tree = 0
  
  for(i in 1:(length(i_seq)-1)) {
    f_tree = i_seq[i]+1
    e_tree = i_seq[i+1]
    
    sampledTreesV_i <- sampledTreesV[seq(f_tree,e_tree,by=tree_every)]
    
    # if(keep_only_ml_events) {
    #
    # rimuovo alberi che includono una back mutations non inclusa nell'albero a ML
    #
    #   re_p_bm <- paste0(possible_bm, collapse = '|')
    #   sampledTreesV_i <- sampledTreesV_i[!grepl(pattern = re_p_bm, x = sampledTreesV_i, fixed = F)]
    #   actual_sampled_tree <- actual_sampled_tree + length(sampledTreesV_i)
    #   if(length(sampledTreesV_i) == 0) {next}
    # }
    
    parentChildPairsV <- unlist(lapply(sampledTreesV_i, function(x){strsplit(x,';')[[1]]}))
    parent_child_pairs <- table(parentChildPairsV)
    
    for(pc_i in 1:length(parent_child_pairs)) {
      parent_child_id <- strsplit(x = names(parent_child_pairs)[pc_i], split = '>', fixed = T)[[1]]
      
      # check if the node is already included in the adjacency matrix otherwise add an empty columns
      if(!(parent_child_id[1] %in% colnames(cons_adjM))) {
        cons_adjM <- cbind(cons_adjM, rep(0, nrow(cons_adjM)))
        cons_adjM <- rbind(cons_adjM, rep(0, ncol(cons_adjM)))
        
        colnames(cons_adjM)[ncol(cons_adjM)] <- rownames(cons_adjM)[nrow(cons_adjM)] <- parent_child_id[1]
      }
      
      if(!(parent_child_id[2] %in% colnames(cons_adjM))) {
        cons_adjM <- cbind(cons_adjM, rep(0, nrow(cons_adjM)))
        cons_adjM <- rbind(cons_adjM, rep(0, ncol(cons_adjM)))
        
        colnames(cons_adjM)[ncol(cons_adjM)] <- rownames(cons_adjM)[nrow(cons_adjM)] <- parent_child_id[2]
      }
      
      cons_adjM[parent_child_id[1],parent_child_id[2]] <- cons_adjM[parent_child_id[1],parent_child_id[2]] + 
        as.numeric(parent_child_pairs[pc_i])
      
    }
    
  }
  
  colnames(cons_adjM)[which(colnames(cons_adjM)=="germline")] <- 'r'
  rownames(cons_adjM)[which(rownames(cons_adjM)=="germline")] <- 'r'
  
  cons_adjM_raw <- cons_adjM
  
  if(keep_only_ml_events) {
    
    # rimuovo colonne con back mutation non incluse nell'albero a ML
    
    # all_bm <-colnames(cons_adjM)[startsWith(x = colnames(cons_adjM), prefix = "-")]
    # 
    # remove_bm <- setdiff(all_bm, ML_back_mutations)
    
    cons_adjM <- cons_adjM[which(!rownames(cons_adjM) %in% possible_bm),
                           which(!colnames(cons_adjM) %in% possible_bm)]
  }
  
  #tot_selected_trees <- actual_sampled_tree
  
  cons_tree <- compute.consensus(cons_adjM = cons_adjM,
                                 tot_tree = tot_selected_trees)
  
  if(!keep_only_ml_events){
    # Prune leaf back mutations observed less than p_rem_bm times
    max_span_adjM <- as.adj.matrix(B = cons_tree$max_span_tree)
    
    leaf_bm <- rownames(max_span_adjM)[which(rowSums(max_span_adjM)==0 & 
                                               startsWith(rownames(max_span_adjM),'-'))]
    
    tree_changed = FALSE
    for(bm in leaf_bm) {
      if(sum(cons_tree$max_span_weigth_adj[,bm]) < p_rem_bm) {
        cons_adjM = cons_adjM[which(rownames(cons_adjM)!=bm),which(colnames(cons_adjM)!=bm)]
        tree_changed = TRUE
      }
    }
    
    if(tree_changed) {
      cons_tree <- compute.consensus(cons_adjM = cons_adjM,
                                     tot_tree = tot_selected_trees)
    }
  }
  
  if(fix_inconsistencies) {
    # Fix inconsistences - back mutations of never happens ones 
    tree_changed = TRUE
    while(tree_changed) {
      
      tree_changed = FALSE
      
      max_span_adjM <- as.adj.matrix(B = cons_tree$max_span_tree)
      bm_IDs <- rownames(max_span_adjM)[startsWith(rownames(max_span_adjM), prefix = "-")]
      
      for(bm in bm_IDs) {
        
        ancestors <- get.ancestor.nodes(adj = max_span_adjM, t_node = bm)[-1]
        
        if(!(bm %in% paste0("-",ancestors))) {
          # back mutation of a never happen ones -> remove
          cons_adjM = cons_adjM[which(rownames(cons_adjM)!=bm),which(colnames(cons_adjM)!=bm)]
          tree_changed = TRUE
        }
        
      }
      
      if(tree_changed) {
        cons_tree <- compute.consensus(cons_adjM = cons_adjM,
                                       tot_tree = tot_selected_trees)
      }
      
    }
    
  }
  
  cons_tree$cons_adjM_raw <- cons_adjM_raw
  
  
  return(cons_tree)
  
}



phylo2clonal <- function(phylo_tree) {
  
  # get the distance between each tip pair
  tips_distances <- ape::cophenetic.phylo(phylo_tree)
  
  nClones = length(phylo_tree$tip.label) 
  #starting from the root (Healthy) connect each tips with the most near one (or more)
  B_adjM <- as.data.frame(matrix(data=0, nrow = nClones, ncol = nClones))
  
  C_tip <- "Healthy"
  
  tips_done = 1
  
  colnames(B_adjM)[tips_done] = C_tip
  rownames(B_adjM)[tips_done] = C_tip
  
  while(tips_done < nClones) {
    
    N_tips = c()
    
    for(Ct in C_tip) {
      
      # determine the most near tips (excluding itself)
      ch_tips <- names(which(tips_distances[Ct,] == min(tips_distances[Ct,][tips_distances[Ct,]>0])))
      ch_tips <- setdiff(ch_tips, colnames(B_adjM))
      
      for(cht in ch_tips) {
        tips_done = tips_done + 1
        colnames(B_adjM)[tips_done] <- cht
        rownames(B_adjM)[tips_done] <- cht
        B_adjM[Ct,cht] <- 1
      }  
      N_tips <- c(N_tips,ch_tips)
    }
    C_tip <- N_tips
  }
  
  colnames(B_adjM) <- gsub(pattern = "Clone", replacement = "", x = colnames(B_adjM))
  colnames(B_adjM)[1] <- 'r'
  rownames(B_adjM) <- colnames(B_adjM)
  
  ord <- c(1,1+order(colnames(B_adjM)[-1]))
  B_adjM <- B_adjM[ord,ord]
  
  B <- as.B.rooted(B_adjM)
  
  return(B)
}
