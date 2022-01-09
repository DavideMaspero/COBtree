FROM rocker/r-ubuntu:20.04

RUN apt-get update -y

RUN apt install -y wget

# Installing R packages
RUN Rscript -e 'install.packages("data.tree",repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("optparse",repos="https://cloud.r-project.org")'

RUN apt install -y openjdk-11-jdk
RUN wget -qO- https://get.nextflow.io | bash \
    && mv nextflow bin/

RUN apt install -y nano

RUN Rscript -e 'install.packages("BiocManager",repos="https://cloud.r-project.org")'
RUN Rscript -e 'BiocManager::install("graph")'
RUN Rscript -e 'BiocManager::install("RBGL")'
RUN Rscript -e 'install.packages("ggplot2",repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("reshape2",repos="https://cloud.r-project.org")'


ADD /COBtree /COBtree/

RUN cd /COBtree/SCITE/ \
&& g++ *.cpp -o scite \
&& cd ../../

RUN cd /COBtree/SCITElik/ \
&& g++ *.cpp -o getSCITEloglik \
&& cd ../../