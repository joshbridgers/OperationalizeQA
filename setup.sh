#!/bin/bash
# Setup bash script for setting up all programs needed for the operationization of QA for Illumina NGS

# Create conda environment and add programs
conda activate OperationalizeQA
conda install -c bioconda bamutil chronqc illumina-interop gatk4 qualimap multiqc picard samtools snpeff

# Install chronQC
git clone https://github.com/nilesh-tawari/ChronQC.git
cd ChronQC
pip install -r requirements.txt
pip install --editable .
