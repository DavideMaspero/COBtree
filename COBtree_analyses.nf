process Generate_topologies {

  publishDir params.basedir + "/Simulation/topologies/perfect_phylogeny", mode: 'link'

  input:
  val n_top from Channel.of( params.n_topologies )
  each n_mut from Channel.fromList( params.mutations )

  output:
  path '*.rds' into sampling_cell_1
  
  """
  Rscript /COBtree/R_scripts/STEP1_topology_generations.R -n $n_top -m $n_mut -o ./ -s $params.seed 
  """
}


process Single_cell_sampling {

  publishDir params.basedir + "/Simulation/datasets/cell_sampled", mode: 'link'

  input:
  path topology from sampling_cell_1.flatMap()
  
  output:
  path '*_branch_[0-9]_sc_*.rds' into pf_genotypes

  """
  Rscript /COBtree/R_scripts/STEP2_single_cell_sampling.R -t $topology -n $params.n_single_cells -r $params.shuffle_cols -p $params.increase_leaf_prob -s $params.seed 
  """
}

process make_SCITE_dataset {
  
  publishDir params.basedir + "/Simulation/datasets/perfect_phylogeny", mode: 'link', pattern: "*.csv"
  
  input:
  tuple val(FP), val(FN), val(M) from Channel.fromList( params.FP_FN_M ) 
  each path(sampled_cells) from pf_genotypes
  
  output:
  tuple path(sampled_cells), path("*.csv") into SCITE_datasets
  path '*.csv' into SCITE_datasets_2
  
  """
  Rscript /COBtree/R_scripts/STEP3_make_input_dataset.R \
          -c $sampled_cells \
          -p $FP -n $FN -m $M \
          -v 3 -t TRUE -q " " \
          -s $params.seed 
  """
  
}

process run_SCITE {

	input:
	tuple path(sampled_cell), path(dataset) from SCITE_datasets
	each mcmc_break from Channel.fromList( params.mcmc_breaks ).flatMap()
 
  output:
  tuple path(sampled_cell), path(dataset), path("*_ml0.gv"), path("*.samples") into SCITE_consensus_tree
 
  script:
 	def n_mut = (dataset =~ "nMut_(.*)_branch")[0][1].toInteger()
	def idx_mut = params.mutations.indexOf(n_mut)
  def mcmc_steps = String.format("%.0f", (mcmc_break.toDouble()*params.mcmc_unit[idx_mut].toDouble()))
  def n_sc = String.format("%.0f", (dataset =~ "sc_(.*)_FP")[0][1].toDouble())
  def FP = String.format("%.5f", (dataset =~ "FP_(.*)_FN")[0][1].toDouble())
  def FN = String.format("%.5f", (dataset =~ "FN_(.*)_M")[0][1].toDouble())
	"""
  /COBtree/SCITE/scite -i $dataset \
                      -n $n_mut \
                      -m $n_sc \
                      -r 10 \
                      -l $mcmc_steps \
                      -fd $FP \
                      -ad $FN \
                      -p 1 \
                      -o ${dataset.baseName+"_MCMC_"+mcmc_steps}
	"""
}

process compute_SCITE_consensus_tree {

  publishDir params.basedir + "/Simulation/consensus_trees", mode: 'link', pattern: "*_consensus.rds" 
  publishDir params.basedir + "/Simulation/metric_files", mode: 'link', pattern: "*_metrics.tsv" 
  
  input:
  tuple path(sampled_cell), path(dataset), path(ml0_file), path(samples_file) from SCITE_consensus_tree

  output:
  path "*_consensus.rds" into collect_SCITE_results
  path "*_metrics.tsv" into collect_SCITE_metrics
  
  """
  Rscript /COBtree/R_scripts/STEP4_compute_SCITE_cons_tree.R -g ${sampled_cell} -d ${dataset} -s ${samples_file} -l ${ml0_file} -p ${params.likelihood_thr}
  """
}

process plot_SCITE_results {

  publishDir params.basedir + "/results", mode: 'copy'
  
  input:
  path metrics from collect_SCITE_metrics.collect()
  
  output:
  path "*.pdf"
  path "metrics.txt"
  
  """
  Rscript /COBtree/R_scripts/STEP5_plot_SCITE_results.R
  """
  
}
