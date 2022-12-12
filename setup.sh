#!/bin/bash

set -eu

###############################################################################
## Setup
###############################################################################

# change to 'hg38' or 'hg19'
REFERENCE_VERSION=hg38
FLOWCELL_ID=MiSeqDemo
FLOWCELL_ID=HW5YNDSX3
DATA_DIR=data
LOG_DIR=logs
RESULTS_DIR=results
LIBRARY_ID=SRR1518133

PerBasePassPhred_PHRED_GTE_CUTOFF=30
OnTargetCovGTPercent_PHRED_GT_CUTOFF=29

#mkdir -v ${DATA_DIR} ${LOG_DIR} ${RESULTS_DIR}

###############################################################################
## Download FASTQ
## Exome of NA12878 (https://www.ncbi.nlm.nih.gov/sra/?term=SRR1518133)
###############################################################################

#wget -o ${LOG_DIR}/fastq_wget.log \
#    -P ${DATA_DIR} \
#    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518133/SRR1518133_1.fastq.gz \
#    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518133/SRR1518133_2.fastq.gz \

###############################################################################
## Download Human Reference  
###############################################################################

## These references are not intended to be used in production and for educational
## use only

## UCSD hg19 genome does NOT have the rCRS version of the mitochondria
hg19_ftp=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

## It is recommend to mask erronous duplications in this reference genome.
## Source: "Failure to Detect Mutations in U2AF1 due to Changes in the GRCh38 
## Reference Sequence", https://doi.org/10.1016/j.jmoldx.2021.10.013
hg38_ftp=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

#wget -o ${LOG_DIR}/reference_genome.log \
#    -P ${DATA_DIR} \
#    $hg19_ftp $hg38_ftp

#wget -o ${LOG_DIR}/snpeff_database_wget.log \
#    -P ${DATA_DIR} \
#    https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_GRCh37.75.zip

###############################################################################
## Download and extract example Illumina SAV files
###############################################################################

## Example from Illumina's InterOp github
## https://github.com/Illumina/interop/blob/master/docs/src/example_sav_analysis.md
#wget -o ${LOG_DIR}/interop_wget.log \
#    -P ${DATA_DIR} \
#    https://github.com/Illumina/interop/releases/download/v1.0.6/MiSeqDemo.zip

#md5sum -c data.md5
#unzip ${DATA_DIR}/MiSeqDemo.zip -d data

###############################################################################
## Install conda
###############################################################################

#wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
#bash Miniconda3-latest-Linux-x86_64.sh -b

###############################################################################
## Create conda environment and add programs
###############################################################################

#conda init bash
## Shell might need to be restarted after init
#conda create --name OperationalizeQA
#conda activate OperationalizeQA
#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda install -c bioconda -y \
#    bamutil fastqc illumina-interop gatk4 hmftools-purple multiqc picard \
#    qualimap samtools snpeff


###############################################################################
## Example aligner and variant caller install
###############################################################################

#conda install -c bioconda -y \
#    bwa minimap2 varscan

###############################################################################
## Install chronQC
###############################################################################

#git clone https://github.com/nilesh-tawari/ChronQC.git
#cd ChronQC
#pip install -r requirements.txt
#pip install --editable .

###############################################################################
## Sequencing Run-Based QC
###############################################################################

## InterOp text output is subject to change with subsequent versions.  Parsing 
## below are meant as examples to be adapted for the users own environment and 
## use case

#interop_summary \
#    ${DATA_DIR}/${FLOWCELL_ID}/ \
#    > ${RESULTS_DIR}/${FLOWCELL_ID}_interop_summary.txt

#interop_plot_qscore_histogram \
#    ${DATA_DIR}/${FLOWCELL_ID} \
#    > ${RESULTS_DIR}/${FLOWCELL_ID}_interop_plot_qscore_histogram.txt

## Cluster density
###############################################################################

## Per Lane cluster density
## head -1 used to get the top surface value.  Depending on the machine the 
## cluster density might be reported as the same across all surfaces. This 
## command will need to be adjusted if one wants to track per lane, per surface
## cluster density
grep "^[0-9]" ${RESULTS_DIR}/${FLOWCELL_ID}_interop_summary.txt |\
    cut -d, -f1,4 |\
    sed 's/ //g' |\
    sort -u |\
    head -1 |\
    awk -F',' -v RESULTS_DIR="$RESULTS_DIR" -v FLOWCELL_ID="$FLOWCELL_ID" '{print $2 > RESULTS_DIR"/"FLOWCELL_ID"_lane_"$1".ClustDen"}'

## Flowcell average cluster density
## dumptext depreciated in recent interop update
#head –n -2 ${RESULTS_DIR}/${FLOWCELL_ID}_interop_dumptext.txt |\
#    sed -n -e '/Lane/,$p' |\
#    tail -n +2 |\
#    awk -F',' '{sum+=$6} END {print sum/NR}' \
#    > "${RESULTS_DIR}/${FLOWCELL_ID}.cluster_density

## Number of reads passing a minimum Phred score criteron
###############################################################################

sed -n '/Level/,/^$/p' ${RESULTS_DIR}/${FLOWCELL_ID}_interop_summary.txt |\
    head -n -1 |\
    awk -F',' '{print $7}' |\
    tail -1 \
    > ${RESULTS_DIR}/${FLOWCELL_ID}.NumReadPassPhred

## Percent of bases higher than the minimum Phred score of all bases called
###############################################################################

## Parses from Illumina's interop summary which uses a minimum Phred score of 
## >=30
egrep ^[0-9]*,[0-9]*\.[0-9]*,1 \
    ${RESULTS_DIR}/${FLOWCELL_ID}_interop_plot_qscore_histogram.txt |\
    awk -F',' -v PHRED_GTE_CUTOFF="$PerBasePassPhred_PHRED_GTE_CUTOFF" '
        {
            if ($1>=PHRED_GTE_CUTOFF) hq_sum+=$2
            else lq_sum+=$2
        } 
            END \
        {
            print hq_sum/(hq_sum + lq_sum) * 100
        }' \
    > ${RESULTS_DIR}/${FLOWCELL_ID}.PerBasePassPhred    

## Demultiplexing success
###############################################################################
grep "Total" ${RESULTS_DIR}/${FLOWCELL_ID}_interop_summary.txt |\
    cut -d, -f3 |\
    sed 's/ \+ /\t/g' \
    > ${RESULTS_DIR}/${FLOWCELL_ID}.DemultiPercent

###############################################################################
## Run Example Alignment
###############################################################################

## For the sake of runtime and brevity steps like BQSR and Indel realignment 
## are skipped

## Build bwa index of the reference
#bwa index ${DATA_DIR}/hg19.fa.gz

## Build index of the reference
#samtools faidx ${DATA_DIR}/hg19.fa.gz

#bwa mem ${DATA_DIR}/hg19.fa.gz \
#    ${DATA_DIR}/${LIBRARY_ID}*_1.fastq.gz \
#    ${DATA_DIR}/${LIBRARY_ID}*_2.fastq.gz \
#    > ${DATA_DIR}/${LIBRARY_ID}.bam

#samtools sort ${DATA_DIR}/${LIBRARY_ID}.sort.bam
#rm ${DATA_DIR}/${LIBRARY_ID}.bam
#samtools index ${DATA_DIR}/${LIBRARY_ID}.sort.bam

#samtools view \
#    -h \
#    -L ${DATA_DIR}/20130108.exome.targets.bed \
#    ${DATA_DIR}/${LIBRARY_ID}.bam \
#    > ${DATA_DIR}/${LIBRARY_ID}.on_target.bam

#picard \
#    MarkDuplicates \
#    -I ${DATA_DIR}/${LIBRARY_ID}.sort.bam \
#    -O ${DATA_DIR}/${LIBRARY_ID}.MarkDuplicates.bam \
#    -M ${DATA_DIR}/${LIBRARY_ID}.MarkDuplicates.txt

###############################################################################
## Run example variant calling and annotation
###############################################################################

#samtools mpileup \
#    -BA \
#    -d 500000 \
#    -q 1 \
#    -f ${DATA_DIR}/hg19.fa.gz \
#    -l ${DATA_DIR}/20130108.exome.targets.bed \
#    ${DATA_DIR}/${LIBRARY_ID}.MarkDuplicates.bam \
#    > ${DATA_DIR}/${LIBRARY_ID}.mpileup

#varscan \
#    mpileup2cns \
#    ${DATA_DIR}/${LIBRARY_ID}.mpileup \
#    --min-var-freq 0.01 \
#    --p-value 0.01 \
#    --strand-filter 0 \
#    --output-vcf \
#    --variants \
#    --min-avg-qual 20 \
#    > ${DATA_DIR}/${LIBRARY_ID}.vcf

#rm ${DATA_DIR}/${LIBRARY_ID}.mpileup

#snpEff \
#    -Xms750m \
#    -csvStats ${DATA_DIR}/${LIBRARY_ID}.effects-stats.csv \
#    GRCh37.75 \
#    ${DATA_DIR}/${LIBRARY_ID}.vcf \
#    > ${DATA_DIR}/${LIBRARY_ID}.annotated.vcf

###############################################################################
## Build interval file 
###############################################################################

#wget -o ${LOG_DIR}/target_bed_wget.log \
#    -P ${DATA_DIR} \
#    ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/exome_pull_down_targets/20130108.exome.targets.bed 

## Create Dictionary file
#picard CreateSequenceDictionary \
#    -R ${DATA_DIR}/hg19.fa.gz

## Create interval file
#picard BedToIntervalList \
#    I=${DATA_DIR}/20130108.exome.targets.bed \
#    SD=${DATA_DIR}/hg19.dict \
#    O=${DATA_DIR}/20130108.exome.targets.interval_list

###############################################################################
## Run Alignment-Based QC
###############################################################################

## Mean on-target coverage of reads
###############################################################################

#picard \
#    CollectHsMetrics \
#    INPUT=${DATA_DIR}/${LIBRARY_ID}.sort.bam \
#    OUTPUT=${DATA_DIR}/${LIBRARY_ID}.CollectHsMetrics.txt \
#    BAIT_INTERVALS=${DATA_DIR}/20130108.exome.targets.interval_list \
#    TARGET_INTERVALS=${DATA_DIR}/20130108.exome.targets.interval_list

grep -A 2 "## METRICS CLASS" \
    ${DATA_DIR}/${LIBRARY_ID}.CollectHsMetrics.txt |\
    cut -f34 |\
    tail -1 \
    > ${RESULTS_DIR}/${LIBRARY_ID}.MeanOnTargetCov

## Percent of targeted bases with coverage greater than a specified minimum
###############################################################################

#gatk DepthOfCoverage \
#    -R ${DATA_DIR}/hg19.fa.gz \
#    -O ${DATA_DIR}/${LIBRARY_ID}.DepthOfCoverage \
#    -I ${DATA_DIR}/${LIBRARY_ID}.sort.bam \
#    -L ${DATA_DIR}/20130108.exome.targets.interval_list \
#    --summary-coverage-threshold ${OnTargetCovGTPercent_PHRED_GT_CUTOFF}

cut -d, -f7 ${DATA_DIR}/${LIBRARY_ID}.DepthOfCoverage.sample_summary \
    > ${RESULTS_DIR}/${LIBRARY_ID}.OnTargetCovGTPercent 

## TODO: Fix depth of coverage crash!

##Percent of bases exceeding the minimum Phred score mapped on target
###############################################################################

#bam stats \
#    --in ${DATA_DIR}/${LIBRARY_ID}.on_target.bam \
#    --phred \
#    2> ${DATA_DIR}/${LIBRARY_ID}.bamutil_stats_on_target.txt

grep ^[0-9] ${DATA_DIR}/${LIBRARY_ID}.bamutil_stats_on_target.txt |\
    awk -v PHRED_GTE_CUTOFF="$PerBasePassPhred_PHRED_GTE_CUTOFF" '
        {
            if ($1>=PHRED_GTE_CUTOFF) hq_sum+=$2
            else lq_sum+=$2
        } 
            END \
        {
            print hq_sum /(hq_sum + lq_sum) * 100
        }'\
    > ${RESULTS_DIR}/${LIBRARY_ID}.OnTargetPercentBase

## AT/GC bias
###############################################################################

#qualimap \
#    bamqc \
#    -bam ${DATA_DIR}/${LIBRARY_ID}.sort.bam

grep "GC percentage" ${DATA_DIR}/${LIBRARY_ID}.sort_stats/genome_results.txt \
    -m1 |\
    tr ' ' '\n' |\
    tail -1 | sed 's/%//' \
    > ${RESULTS_DIR}/${LIBRARY_ID}.GCBias

## Median insert size (bp)
###############################################################################

#picard CollectInsertSizeMetrics \
#    -I ${DATA_DIR}/${LIBRARY_ID}.sort.bam \
#    -O ${DATA_DIR}/${LIBRARY_ID}.CollectInsertSizeMetrics.txt \
#    -H ${DATA_DIR}/${LIBRARY_ID}.CollectInsertSizeMetrics.pdf \

grep -A2 '## METRICS CLASS' \
    ${DATA_DIR}/${LIBRARY_ID}.CollectInsertSizeMetrics.txt |\
    cut -f1 |\
    tail -1 \
    > ${RESULTS_DIR}/${LIBRARY_ID}.MeanInsertSize

## Percent duplicates
###############################################################################

grep -A2  "## METRICS CLASS" \
    ${DATA_DIR}/${LIBRARY_ID}.MarkDuplicates.txt |\
    cut -f9 |\
    tail -1 \
    > ${RESULTS_DIR}/${LIBRARY_ID}.PercentDuplicate

## Observed sex matches reported sex
###############################################################################

#egrep "[XY]:" <SAMPLE_NAME>.DepthOfCoverage | sed 's/:/,/' | awk -F',' '{if ($1 ~ /X$/) {X_sum+=4 ; X_count+=1} else {Y_sum+=4 ; Y_count+=1}} END {print "("X_sum"/"X_count")/("Y_sum"/"Y_count")"}' | bc -l

#    > ${RESULTS_DIR}/${LIBRARY_ID}.SexMatch

###############################################################################
# Run Variant-Based QC
###############################################################################

## SNV/INDEL ratio
###############################################################################

grep "^SNP ," ${DATA_DIR}/${LIBRARY_ID}.effects-stats.csv |\
    cut -d, -f3 |\
    sed 's/[ %]//g' \
    > ${RESULTS_DIR}/${LIBRARY_ID}.SNVINDELRatio

## Ti/TV ratio
###############################################################################

grep "^Ts_Tv_ratio" ${DATA_DIR}/${LIBRARY_ID}.effects-stats.csv  |\
    cut -d, -f2 |\ 
    sed 's/[ %]//g' \ 
    > ${RESULTS_DIR}/${LIBRARY_ID}.TiTvRatio
