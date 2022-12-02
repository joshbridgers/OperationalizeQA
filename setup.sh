#!/bin/bash

set -eu

###############################################################################
# Setup
###############################################################################

# change to 'hg38' or 'hg19'
REFERENCE_VERSION=hg38

#mkdir -v data logs results 

###############################################################################
# Download FASTQ
# Exome of NA12878 (https://www.ncbi.nlm.nih.gov/sra/?term=SRR1518133)
###############################################################################

#wget -o logs/fastq_wget.log \
#    -P data \
#    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518133/SRR1518133_1.fastq.gz \
#    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518133/SRR1518133_2.fastq.gz 

###############################################################################
# Download Human Reference  
###############################################################################

# These references are not intended to be used in production and for educational
# use only

hg19_ftp=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
# UCSD hg19 genome does NOT have the rCRS version of the mitochondria
hg38_ftp=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
# It is recommend to mask erronous duplications in this reference genome.
# Source: "Failure to Detect Mutations in U2AF1 due to Changes in the GRCh38 
# Reference Sequence", https://doi.org/10.1016/j.jmoldx.2021.10.013
#wget -o logs/reference_genome.log \
#    -P data \
#    $hg19_ftp $hg38_ftp

###############################################################################
# Download and extract example Illumina SAV files
###############################################################################

# Example from Illumina's InterOp github
# https://github.com/Illumina/interop/blob/master/docs/src/example_sav_analysis.md
#wget -o logs/interop_wget.log \
#    -P data \
#    https://github.com/Illumina/interop/releases/download/v1.0.6/MiSeqDemo.zip

#md5sum -c data.md5
#unzip data/MiSeqDemo.zip -d data

###############################################################################
#Install conda
###############################################################################

#wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
#bash Miniconda3-latest-Linux-x86_64.sh -b

###############################################################################
#Create conda environment and add programs
###############################################################################

#conda init bash
# Shell might need to be restarted after init
#conda create --name OperationalizeQA
#conda activate OperationalizeQA
#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda install -c bioconda -y \
#    bamutil fastqc illumina-interop gatk4 hmftools-purple multiqc picard \
#    qualimap samtools snpeff


###############################################################################
#Aligner install
###############################################################################

#conda install -c bioconda -y \
#    bwa

###############################################################################
# Install chronQC
###############################################################################
#git clone https://github.com/nilesh-tawari/ChronQC.git
#cd ChronQC
#pip install -r requirements.txt
#pip install --editable .

###############################################################################
# Sequencing Run QC
###############################################################################

# Cluster density
interop_summary \
    data/MiSeqDemo/ \
    > example.interop_summary.txt

grep "Read [0-9]$" -A2 example.interop_summary.txt | grep ",Density" -A1 | grep ^[0-9] | cut -d, -f1,4 | sed 's/ //g' | awk -F',' '{print $2 > "results/Lane_"$1"_cluster_density"}'

###############################################################################
# Run Alignment
###############################################################################

###############################################################################
# Run Alignment QC
###############################################################################
