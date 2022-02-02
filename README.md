# COBtree
Improve the convergence of MCMC-based inference framework  of mutation trees applying the Optimum Branching algorithm on the trees sampled from the posterior


## Installation
1. Clone the COBtree repository
```
git clone https://github.com/DavideMaspero/COBtree.git
```
2. Install Docker following the instruction at this link: (https://www.docker.com/get-started).
3. Get COBtree docker image. There are two options to do this:
   - Download it from Docker Hub (https://hub.docker.com/u/dcblab)
     ```
     docker pull dcblab/cobtree_img:latest
     ```
   - Create it from docker file
     ```
     cd docker
     docker build -t cobtree_img .
     ```
## Analysis
0. Parameters can be changed by editing the *nextflow.config* file in COBtree directory (plese, do not rename it)
1. Run the docker image iteratively by specifing as source directory the local **absolute** path where COBtree repository is located. The SCITE inferences are executed in parallel. So, please replace *n_cpus* with the number of CPUs available. 
```
docker run -it --rm --cpus n_cpus --mount type=bind,source=/local-dir-absolute-path/COBtree/,target=/nextflow_runs/ cobtree_img:latest
```
2. Inside the cobtree_img, move to the nextflow_run directory
```
cd nextflow_run
```
3. Run the Simulation
```
nextflow run COBtree_analyses.nf
```
## Results
The results (plots and metric scores) are stored in the COBtree/results directory.  
