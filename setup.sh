#!/bin/bash
# Example implementation of parsing for key QC metrics from common 
# bioinformatics pipeline on human cancer samples 
#
################################################################################
# WARNING
################################################################################
# LIMITATIONS
#
# Commands are parsing text output and can break if the format changes. 
# Commands have not been tuned for multiple cores and operated under default
# behaviors. These scripts are NOT directly intended for clinical use. All 
# production code will need to be validated for use in a clinical setting. 

# "strict" mode
# https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425
set -euxo pipefail

################################################################################
# Setup
################################################################################
DATA_DIR=data
LOG_DIR=logs
GATK3_8_LOC=https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
MINICONDA_LOC=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Functions
function wget_func () { 
    wget \
        -P ${DATA_DIR} \
        $1 2>&1 | tee -a ${LOG_DIR}/${2}_wget.log
}

mkdir -vp ${DATA_DIR} ${LOG_DIR}

################################################################################
# Setup   
################################################################################

# Download conda
################################################################################
wget_func ${MINICONDA_LOC} "miniconda"

# Perform checksum on downloads
################################################################################
md5sum -c setup_files.md5

# Install conda
################################################################################
bash ${DATA_DIR}/Miniconda3-latest-Linux-x86_64.sh -b

# Create conda environment and add programs
################################################################################
conda create --name OperationalizeQA
conda init bash
# Shell might need to be restarted after init
conda activate OperationalizeQA
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda -y \
    bamutil fastqc illumina-interop gatk4 hmftools-purple multiqc picard \
    qualimap samtools snpeff

# Example aligners and variant caller install
conda install -c bioconda -y \
    bwa bwa-mem2 minimap2 varscan

# Alternative GATK3 install
################################################################################
#wget_func ${GATK3_8_LOC} "gatk3.8"
#GATK3_8_LOC=$(basename ${GATK3_8_LOC})
#bunzip ${GATK_NAME} 
#tar xfv ${GATK_NAME} 

# Install chronQC
################################################################################
git clone https://github.com/nilesh-tawari/ChronQC.git
cd ChronQC
pip install -r requirements.txt
pip install --editable .
