#!/usr/bin/env bash

# Example implementation of parsing for key QC metrics from common 
# bioinformatics pipeline on human cancer samples 
#
################################################################################
# Warnings and Limitations
################################################################################
# See https://github.com/joshbridgers/OperationalizeQA#warnings-and-limitations

################################################################################
# Setup
################################################################################

# "strict" mode
set -euxo pipefail

# Directory locations
DATA_DIR=data
FLOWCELL_ID=MiSeqDemo 
LIBRARY_ID=SRR1518133
LOG_DIR=logs
RESULTS_DIR=results

# Thresholds
PerBasePassPhred_PHRED_GTE_CUTOFF=30
OnTargetCovGTPercent_PHRED_GT_CUTOFF=29

# Download locations

# Publically available human cell line NA12878
FASTQ_1_LOC=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518133/SRR1518133_1.fastq.gz
FASTQ_2_LOC=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518133/SRR1518133_2.fastq.gz
INTEROP_LOC=https://github.com/Illumina/interop/releases/download/v1.0.6/MiSeqDemo.zip
REFERENCE_LOC=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
REFERENCE_INDEX_LOC=https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai
REFERENCE_FILE_SUFFIX=$(basename ${REFERENCE_LOC} | sed 's/\.fasta\.gz$//')
SNPEFF_DATABASE_LOC=https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_GRCh37.75.zip
TARGET_BED_LOC=ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/exome_pull_down_targets/20130108.exome.targets.bed 
TARGET_BED_SUFFIX=$(basename ${TARGET_BED_LOC} | sed 's/\.bed//')

mkdir -vp ${DATA_DIR} ${LOG_DIR} ${RESULTS_DIR}

################################################################################
# Functions
################################################################################
function wget_func () { 
    wget \
        -P ${DATA_DIR} \
        $1 2>&1 | tee -a ${LOG_DIR}/${2}_wget.log
}

################################################################################
# Downloads   
################################################################################

# Download human reference and index
################################################################################
# GRCh37 
# A no-alt version was used. This reference does not include decoy sequence
# hs37d5.  

# hg38
# It is recommend to mask erronous duplications in this reference genome.
# Source: "Failure to Detect Mutations in U2AF1 due to Changes in the GRCh38 
#   Reference Sequence", https://doi.org/10.1016/j.jmoldx.2021.10.013. Accessed 
#   2022-12-12

wget_func ${REFERENCE_LOC} "reference"
wget_func ${REFERENCE_INDEX_LOC} "reference_index"

# Download FASTQ
################################################################################
# Exome of NA12878 (https://www.ncbi.nlm.nih.gov/sra/?term=SRR1518133)
wget_func $FASTQ_1_LOC "fastq_1"
wget_func $FASTQ_2_LOC "fastq_2"

# Download snpeff annotation database
################################################################################
wget_func ${SNPEFF_DATABASE_LOC} "snpeff_database"

# Downloads and extract example Illumina SAV files
################################################################################
# Example from Illumina's InterOp 
# https://github.com/Illumina/interop/blob/master/docs/src/example_sav_analysis.md
# This interOp folder is already present in the git repo under the "data" folder
#wget_func ${INTEROP_LOC} "interop"
#unzip ${DATA_DIR}/MiSeqDemo.zip -d data

# Download bed used by example exome
################################################################################
wget_func ${TARGET_BED_LOC} "exome_bed"

# Perform checksum on downloads
################################################################################
md5sum -c ${DATA_DIR}/run_files.md5

# Remove "chr" from bed
################################################################################
perl -pi -e 's/^chr//' ${DATA_DIR}/${TARGET_BED_SUFFIX}.bed

# Trim Reference Genome
################################################################################
# truncating needed due to the following error from samtools for this specific
# fastq.gz
# [E::bgzf_read_init] Cannot decompress legacy RAZF format
# [E::razf_info] To decompress this file, use the following commands:
#    truncate -s 891946027 data/human_g1k_v37.fasta.gz
#        gunzip data/human_g1k_v37.fasta.gz
#        The resulting uncompressed file should be 3153506519 bytes in length.
#        If you do not have a truncate command, skip that step (though gunzip will
#        likely produce a "trailing garbage ignored" message, which can be ignored).

truncate -s 891946027 ${DATA_DIR}/${REFERENCE_FILE_SUFFIX}.fasta.gz
gunzip ${DATA_DIR}/${REFERENCE_FILE_SUFFIX}.fasta.gz

################################################################################
# Sequencing Run-Based QC
################################################################################

# InterOp text output is subject to change with subsequent versions.  Parsing 
# below are meant as examples to be adapted for the users own environment and 
# use case

interop_summary \
    ${DATA_DIR}/${FLOWCELL_ID} \
    > ${DATA_DIR}/${FLOWCELL_ID}_interop_summary.txt

interop_plot_qscore_histogram \
    ${DATA_DIR}/${FLOWCELL_ID} \
    > ${DATA_DIR}/${FLOWCELL_ID}_interop_plot_qscore_histogram.txt

# Cluster density
################################################################################
# Per Lane cluster density
# head -1 used to get the top surface value.  Depending on the machine the 
# cluster density might be reported as the same across all surfaces. This 
# command will need to be adjusted if one wants to track per lane, per surface
# cluster density. Illumina reports the average cluster density and the 
# standard deviations across all tiles for that lane.
grep "^[0-9]" ${DATA_DIR}/${FLOWCELL_ID}_interop_summary.txt |\
    cut -d, -f1,4 |\
    sed 's/ //g' |\
    sort -u |\
    head -1 |\
    awk -F',' -v RESULTS_DIR="$RESULTS_DIR" \
        -v FLOWCELL_ID="$FLOWCELL_ID" '
        {
            print $2 > RESULTS_DIR"/"FLOWCELL_ID"_lane_"$1".ClustDen"
        }'

# Flowcell average cluster density
# Use the awk command below for an average across the flowcell  
#    awk -F',' '{sum+=$2} END {print sum/NR}' \
#    > ${RESULTS_DIR}/${FLOWCELL_ID}.ClustDen

# Number of reads passing a minimum Phred score criteron
################################################################################
sed -n '/Level/,/^$/p' ${DATA_DIR}/${FLOWCELL_ID}_interop_summary.txt |\
    head -n -1 |\
    awk -F',' '{print $7}' |\
    tail -1 \
    > ${RESULTS_DIR}/${FLOWCELL_ID}.NumReadPassPhred

# Percent of bases higher than the minimum Phred score of all bases called
################################################################################
# Parses from Illumina's interop summary which uses a minimum Phred score of 
# >=30
egrep ^[0-9]*,[0-9]*\.[0-9]*,1 \
    ${DATA_DIR}/${FLOWCELL_ID}_interop_plot_qscore_histogram.txt |\
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

# Demultiplexing success
################################################################################
grep "Total" ${DATA_DIR}/${FLOWCELL_ID}_interop_summary.txt |\
    cut -d, -f3 |\
    sed 's/ \+ /\t/g' \
    > ${RESULTS_DIR}/${FLOWCELL_ID}.DemultiPercent

################################################################################
# Run Example Alignment
################################################################################
# For the sake of runtime and brevity steps like BQSR and Indel realignment 
# are skipped

bwa index ${DATA_DIR}/${REFERENCE_FILE_SUFFIX}.fasta
bwa mem \
    ${DATA_DIR}/${REFERENCE_FILE_SUFFIX}.fasta \
    -R "@RG\tID:${LIBRARY_ID}\tPL:illumina\tPU:${LIBRARY_ID}\tSM:${LIBRARY_ID}" \
    ${DATA_DIR}/${LIBRARY_ID}*_1.fastq.gz \
    ${DATA_DIR}/${LIBRARY_ID}*_2.fastq.gz \
    > ${DATA_DIR}/${LIBRARY_ID}.bam

samtools sort ${DATA_DIR}/${LIBRARY_ID}.bam \
    -o ${DATA_DIR}/${LIBRARY_ID}.sort.bam
rm ${DATA_DIR}/${LIBRARY_ID}.bam

samtools view \
    -h \
    -L ${DATA_DIR}/${TARGET_BED_SUFFIX}.bed \
    ${DATA_DIR}/${LIBRARY_ID}.sort.bam \
    > ${DATA_DIR}/${LIBRARY_ID}.on_target.bam
rm ${DATA_DIR}/${LIBRARY_ID}.sort.bam

picard \
    MarkDuplicates \
    -I ${DATA_DIR}/${LIBRARY_ID}.on_target.bam \
    -O ${DATA_DIR}/${LIBRARY_ID}.MarkDuplicates.bam \
    -M ${DATA_DIR}/${LIBRARY_ID}.MarkDuplicates.txt
samtools index ${DATA_DIR}/${LIBRARY_ID}.MarkDuplicates.bam

################################################################################
# Run example variant calling and annotation
################################################################################
samtools mpileup \
    -BA \
    -d 500000 \
    -q 1 \
    -f ${DATA_DIR}/${REFERENCE_FILE_SUFFIX}.fasta \
    -l ${DATA_DIR}/${TARGET_BED_SUFFIX}.bed \
    ${DATA_DIR}/${LIBRARY_ID}.MarkDuplicates.bam \
    > ${DATA_DIR}/${LIBRARY_ID}.mpileup

varscan \
    mpileup2cns \
    ${DATA_DIR}/${LIBRARY_ID}.mpileup \
    --min-var-freq 0.01 \
    --p-value 0.05 \
    --strand-filter 0 \
    --output-vcf \
    --variants \
    --min-avg-qual 20 \
    > ${DATA_DIR}/${LIBRARY_ID}.vcf

rm ${DATA_DIR}/${LIBRARY_ID}.mpileup

snpEff \
    -Xms750m \
    -csvStats ${DATA_DIR}/${LIBRARY_ID}.effects-stats.csv \
    -stats ${DATA_DIR}/${LIBRARY_ID}_snpEff_summary.html \
    GRCh37.75 \
    ${DATA_DIR}/${LIBRARY_ID}.vcf \
    > ${DATA_DIR}/${LIBRARY_ID}.annotated.vcf

################################################################################
# Build dictionary and interval file 
################################################################################
picard CreateSequenceDictionary \
    -R ${DATA_DIR}/${REFERENCE_FILE_SUFFIX}.fasta

picard BedToIntervalList \
    I=${DATA_DIR}/${TARGET_BED_SUFFIX}.bed \
    SD=${DATA_DIR}/${REFERENCE_FILE_SUFFIX}.dict \
    O=${DATA_DIR}/${TARGET_BED_SUFFIX}.interval_list

################################################################################
# Run alignment-Based QC
################################################################################

# Mean on-target coverage of reads
################################################################################
# Note a probe bed was not provided in this example so the target bed was reused
# for the BAIT_INTERVALS
picard \
    CollectHsMetrics \
    INPUT=${DATA_DIR}/${LIBRARY_ID}.MarkDuplicates.bam \
    OUTPUT=${DATA_DIR}/${LIBRARY_ID}.CollectHsMetrics.txt \
    BAIT_INTERVALS=${DATA_DIR}/${TARGET_BED_SUFFIX}.interval_list \
    TARGET_INTERVALS=${DATA_DIR}/${TARGET_BED_SUFFIX}.interval_list

grep -A 2 "## METRICS CLASS" \
    ${DATA_DIR}/${LIBRARY_ID}.CollectHsMetrics.txt |\
    cut -f34 |\
    tail -1 \
    > ${RESULTS_DIR}/${LIBRARY_ID}.MeanOnTargetCov

# Percent of targeted bases with coverage greater than a specified minimum
################################################################################
gatk DepthOfCoverage \
    -R ${DATA_DIR}/${REFERENCE_FILE_SUFFIX}.fasta \
    -O ${DATA_DIR}/${LIBRARY_ID}.DepthOfCoverage \
    -I ${DATA_DIR}/${LIBRARY_ID}.MarkDuplicates.bam \
    -L ${DATA_DIR}/${TARGET_BED_SUFFIX}.interval_list \
    --summary-coverage-threshold ${OnTargetCovGTPercent_PHRED_GT_CUTOFF}

grep ${LIBRARY_ID} ${DATA_DIR}/${LIBRARY_ID}.DepthOfCoverage.sample_summary |\
    cut -d, -f7 \
    > ${RESULTS_DIR}/${LIBRARY_ID}.OnTargetCovGTPercent 

#Percent of bases exceeding the minimum Phred score mapped on target
################################################################################
bam stats \
    --in ${DATA_DIR}/${LIBRARY_ID}.on_target.bam \
    --phred \
    2> ${DATA_DIR}/${LIBRARY_ID}.bamutil_stats_on_target.txt

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

# AT/GC bias
################################################################################
qualimap \
    bamqc \
    -bam ${DATA_DIR}/${LIBRARY_ID}.MarkDuplicates.bam

grep "GC percentage" ${DATA_DIR}/${LIBRARY_ID}.MarkDuplicates_stats/genome_results.txt \
    -m1 |\
    tr ' ' '\n' |\
    tail -1 | sed 's/%//' \
    > ${RESULTS_DIR}/${LIBRARY_ID}.GCBias

# Median insert size (bp)
################################################################################
picard CollectInsertSizeMetrics \
    -I ${DATA_DIR}/${LIBRARY_ID}.MarkDuplicates.bam \
    -O ${DATA_DIR}/${LIBRARY_ID}.CollectInsertSizeMetrics.txt \
    -H ${DATA_DIR}/${LIBRARY_ID}.CollectInsertSizeMetrics.pdf \

grep -A2 '## METRICS CLASS' \
    ${DATA_DIR}/${LIBRARY_ID}.CollectInsertSizeMetrics.txt |\
    cut -f1 |\
    tail -1 \
    > ${RESULTS_DIR}/${LIBRARY_ID}.MeanInsertSize

# Percent duplicates
################################################################################
grep -A2  "## METRICS CLASS" \
    ${DATA_DIR}/${LIBRARY_ID}.MarkDuplicates.txt |\
    cut -f9 |\
    tail -1 \
    > ${RESULTS_DIR}/${LIBRARY_ID}.PercentDuplicate

# Observed sex matches reported sex
################################################################################
egrep "[XY]:" ${DATA_DIR}/${LIBRARY_ID}.DepthOfCoverage |\
    sed 's/:/,/' |\
    awk -F',' '
        {
            if ($1 ~ /X$/) {
                X_sum+=4 ; X_count+=1
            } 
            else {
                Y_sum+=4 ; Y_count+=1
            }
        } 
            END \
        {
            print "("X_sum"/"X_count")/("Y_sum"/"Y_count")"
        }' |\
    bc -l \
    > ${RESULTS_DIR}/${LIBRARY_ID}.SexMatch

################################################################################
# Run variant-based QC
################################################################################

# SNV/INDEL ratio
################################################################################
grep "^SNP ," ${DATA_DIR}/${LIBRARY_ID}.effects-stats.csv |\
    cut -d, -f3 |\
    sed 's/[ %]//g' \
    > ${RESULTS_DIR}/${LIBRARY_ID}.SNVINDELRatio

# Ti/TV ratio
################################################################################
grep "^Ts_Tv_ratio" ${DATA_DIR}/${LIBRARY_ID}.effects-stats.csv  |\
    cut -d, -f2 |\
    sed 's/[ %]//g' \
    > ${RESULTS_DIR}/${LIBRARY_ID}.TiTvRatio
