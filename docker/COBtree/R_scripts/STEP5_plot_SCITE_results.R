############## Read input parameters
suppressPackageStartupMessages(library("optparse"))

option_list <- list(              
  make_option(c("-i", "--input_dir"), type="character", default = "./", 
              help="Directory where TSV metric files, computed by STEP5, are stored",
              metavar="string"),
  
  make_option(c("-o", "--output_dir"), type="character", default = "./",
              help="Path of the directory in which save the results.",
              metavar="string")
)

opt <- parse_args(OptionParser(option_list=option_list))

tsv_fn_list <- dir(path=opt$input_dir, pattern = "_metrics.tsv", include.dirs = F)

if(sum(file.exists(tsv_fn_list))<2) {
  stop("Less than two files detected, please check input file list")
}

out_dir <- opt$output_dir
if(!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = T, showWarnings = F)
}

library(ggplot2)

for(i in 1:length(tsv_fn_list)) {
  tsv_i <- read.table(file = tsv_fn_list[i], header=T, sep='\t')
  if(i == 1) {
    metrics <- tsv_i
  }  
  
  metrics <- rbind(metrics, tsv_i)
}

metrics$noise <- 'Low'
metrics$noise[metrics$missing == 0.1] <- 'Middle'
metrics$noise[metrics$missing == 0.2] <- 'High'
metrics$noise <- factor(x = metrics$noise, levels = c('Low', 'Middle','High'))

metrics$MCMClength <- NA

n_muts <- unique(metrics$n_mut)

for(nMs_i in n_muts) {
  MCMCsteps_us <- sort(unique(metrics$MCMCsteps[metrics$n_mut == nMs_i]))
  
  for(i in 1:length(MCMCsteps_us)){
    metrics$MCMClength[metrics$n_mut == nMs_i & metrics$MCMCsteps == MCMCsteps_us[i]] <- i
  }
  
}

lvl <- sort(as.numeric(unique(metrics$MCMClength)))

metrics$MCMClength <- factor(metrics$MCMClength, levels = lvl, labels = paste0('MCMC len ', lvl))

metrics$deltaPC <- metrics$ML_PCdist - metrics$COB_PCdist
metrics$deltaCG <- metrics$ML_CGerrors - metrics$COB_CGerrors

for(lthr in unique(metrics$lik_thr)) {

pdf(paste0(out_dir,"/CG_PC_scatter_COBthr_",lthr,".pdf"), width = 20, height = 10)
print(
ggplot(data = metrics[metrics$lik_thr==lthr,], mapping = aes(x = deltaPC, y = deltaCG, color = as.factor(noise)))+
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) +
  geom_point()+
  facet_grid(n_mut ~ MCMClength, scale = "free")+
  theme_bw()
)
dev.off()

reshaped_PC <- reshape2::melt(metrics[metrics$lik_thr==lthr, c("SCITEfn",
                                         "n_mut",
                                         "MCMClength",
                                         "noise",
                                         "ML_PCdist",
                                         "COB_PCdist")],
                             id.vars = c("SCITEfn", "MCMClength","n_mut", "noise"),
                             value.name = "Errors")
colnames(reshaped_PC) <- c("Exp_name","MCMClength","n_mut","noise","Tree_type","PC_dist")

reshaped_PC$Tree_type <- factor(x = as.character(reshaped_PC$Tree_type), 
                                levels = c("ML_PCdist", "COB_PCdist"),
                                labels = c("MLtree", "COBtree")) 
                                
                                
                                
reshaped_CG <- reshape2::melt(metrics[metrics$lik_thr==lthr, c("SCITEfn",
                                          "n_mut",
                                          "MCMClength",
                                          "noise",
                                          "ML_CGerrors",
                                          "COB_CGerrors")],
                              id.vars = c("SCITEfn", "MCMClength","n_mut", "noise"),
                              value.name = "Errors")
colnames(reshaped_CG) <- c("Exp_name","MCMClength","n_mut","noise","Tree_type","CG_errors")

reshaped_CG$Tree_type <- factor(x = as.character(reshaped_CG$Tree_type), 
                                levels = c("ML_CGerrors", "COB_CGerrors"),
                                labels = c("MLtree", "COBtree")) 




pdf(paste0(out_dir,"/scite_ML_vs_COB_PCdist_distribution_COBthr_",lthr,".pdf"), width = 25, height = 15)
print(
ggplot(data = reshaped_PC, mapping = aes(x = Tree_type, y = PC_dist, fill=Tree_type)) + 
  geom_violin(draw_quantiles = 0.5) + 
  facet_grid(n_mut ~ MCMClength, scales = "free_y", labeller = label_both)+
  labs(fill="SCITE tree")+
  xlab("")+
  theme_bw()+
  theme(text = element_text(size = 25)) 
)
dev.off()


pdf(paste0(out_dir,"/scite_ML_vs_COB_CGerror_distribution_COBthr_",lthr,".pdf"), width = 25, height = 15)
print(
ggplot(data = reshaped_CG, mapping = aes(x = Tree_type, y = CG_errors, fill=Tree_type)) + 
  geom_violin(draw_quantiles = 0.5) + 
  facet_grid(n_mut ~ MCMClength, scales = "free_y", labeller = label_both)+
  labs(fill="SCITE tree")+
  xlab("")+
  theme_bw()+
  theme(text = element_text(size = 25)) 
)
dev.off()

}


write.table(file=paste0(out_dir,"/metrics.txt"),x=metrics, sep='\t',row.names=F, col.names=T, quote=F)