source("/COBtree/R_functions/utils.R")

############## Read input parameters
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
   make_option(c("-n", "--topologies"), type="integer", default=20, 
               help="Number of random topologies which will be generated [default %default]",
               metavar="number"),
   make_option(c("-m", "--mutations"), type="integer", default=50, 
               help="Number of mutations (nodes) included in each topology [default %default]",
               metavar="number"),
   make_option(c("-o", "--output_dir"), type="character", default = "topologies",
               help="Path of the output directory in which the topologies will be saved [default %default]",
               metavar="string"),
   make_option(c("-b", "--branch_min"), type="integer", default = 2,
               help="Lower bound of the branching size range. [default %default]",
               metavar="number"),
   make_option(c("-B", "--branch_max"), type="integer", default = 5,
               help="Upper bound of the branching size range. [default %default]",
               metavar="number"),
   make_option(c("-s", "--seed"), type="integer", default = sample(1:10000, 1),
               help="Seed for random number, is used for reproducibility. [default: random]",
               metavar="number")
   
)

opt <- parse_args(OptionParser(option_list=option_list))

outdir <- normalizePath(opt$output_dir, winslash = '\\')

nMut <- opt$mutations

nRun <- opt$topologies

b_min <- opt$branch_min
b_max <- opt$branch_max

set.seed(opt$seed)

##############

n_top_offset <- 0

if(!dir.exists(outdir)) {
  dir.create(path = outdir, recursive = T)
} else {
   n_top_offset <- length(dir(path = outdir, pattern = "topology_"))
}

if(nRun > 1) {
   branching_v <- sort(floor(linMap(x = rnorm(nRun), from = b_min, to = b_max)))
   
} else {
   branching_v <- round((b_min + b_max) / 2)
}

for(i in 1:nRun) {
  # randomly set the max brach size (i.e., number of node child)
  range_children = 1:branching_v[i]
  
  # Generat ground-truth tree
  Bgt = matrix(0, nrow=(nMut+1), ncol = (nMut+1))
  diag(Bgt) <- 1
  
  node_left <- 2:(nMut+1)
  
  for(r in 1:nMut) {
    
    if(length(range_children)>1) {
      
      n_child = sample(range_children, size=1, prob=rev(range_children)/sum(range_children))
      
    } else {
      
      n_child = range_children
    }
    
    nc = 0
    
    while(length(node_left)>0 & nc < n_child) {
      nc = nc + 1
      Bgt[node_left[1], 1:r] <- Bgt[r,1:r]
      node_left <- node_left[-1]
    }
  }
  
  #Bgt <- cbind(rep(1,nrow(unroot_Bgt)), unroot_Bgt)
  #Bgt <- rbind(c(1,rep(0,(ncol(Bgt)-1))), Bgt)
  colnames(Bgt) <- c('r',1:nMut)
  rownames(Bgt) <- colnames(Bgt)
  
  topology <- list(Bgt = Bgt, clonal_genotypes = Bgt[-1,-1], branch = branching_v[i])
  
  fn <- paste0(outdir,"/topology_",i+n_top_offset,"_nMut_",nMut,"_branch_",branching_v[i],".rds")
  
  saveRDS(topology, file = fn)
}

