# Operationalizing Quality Assurance for Clinical Illumina Somatic NGS Pipelines

Example implementation of parsing for key QC metrics from a common 
bioinformatics pipeline on human cancer samples.

## Quick Start
##### Setup programs and `conda` environment
```./setup.sh```
##### Run example pipeline 
```./run.sh```

## Documentation
  
### Table 1. Summary of run-level metrics with related bioinformatics software and modules for generation

|Metric|Tool and Module|
|---|---|
|Cluster density|**InterOp** summary|
|Number of reads passing a minimum Phred score criterion|**InterOp** summary|
|Percent of bases higher than the minimum Phred score of all bases called|**InterOp** plot_qscore_histogram|
|Demultiplexing success|**InterOp** index-summary|

### Table 2. Summary of per-sample metrics with the related bioinformatics software and modules for generation

|Metric|Tool and Module|bcbio|
|---|---|---|
|Mean on-target coverage of reads|**Picard** CollectHsMetrics|Yes|
|Percent of targeted bases with coverage greater than a specified minimum|**GATK-4.X** DepthOfCoverage|No|
|Percent of bases exceeding the minimum Phred score mapped on target|**BamUtil** stats|No|
|AT/GC bias|**Qualimap** bamqc|Yes|
|Mean insert size (bp)|**Samtools** stats|Yes|
|Percent PCR duplicates|**Picard** MarkDuplicates|Yes|
|Observed sex matches reported sex|**GATK-4.X** DepthOfCoverage|Yes|
|SNV/INDEL ratio|**SnpEff** csvStats|Yes|
|Ti/Tv ratio|**SnpEff** csvStats|Yes|

## Warnings and Limitations
All data used is publically available.

Commands have not been optimal for multiple cores or balanced in a high 
performance computing environment. Programs will operated under default 
behaviors which may be sub-optimal for said environments. 

These scripts and genome reference are NOT directly intended for clinical 
use but for educational use only. All production code will need to be 
validated for use in a clinical setting. 

The provided commands are parsing text output and can break if the outputted 
format changes. For example, given the constant evolution of the Illumina 
platform, there have been different outputs from the `interop` program that 
range from depreciated modules, to formatting with nuanced differences such 
as per flowcell surface cluster densities as highlighted below.

Although the `conda` command will grab the latest version of the inputted 
program, this example covers parsing for InterOp v1.2.0 release and scripts 
may have to be adapted for older versions of the InterOp folder or program 
output.

All testing was performed on CentOS Linux 7 (Core) with Kernel Linux 
3.10.0-957.5.1.el7.x86_64.
