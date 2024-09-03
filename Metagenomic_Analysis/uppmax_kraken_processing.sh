#!/bin/bash -l

# tell it is bash language and -l is for starting a session with a "clean environment, e.g. with no modules loaded and paths reset"

#SBATCH -A naiss2023-22-1288  # Project name

#SBATCH -p core  # Asking for cores (for test jobs and as opposed to multiple nodes) 

#SBATCH -n 16  # Number of cores

#SBATCH -t 24:00:00  #

#SBATCH -J M14-samples-kraken-results-processing  # Name of the job

# go to some directory

cd /proj/naiss2023-23-618/FIBER_INT_METAGENOMIC_DATA/02_RAW_READS/CLEAN_READS_V1/TEST_KRAKEN_PROCESSING/M159_DC
pwd -P

# load software modules

module load bioinfo-tools
module load bioinfo-tools conda
module load bioinfo-tools bwa/0.7.18
module load bioinfo-tools megahit/1.2.9
module load bioinfo-tools SeqKit/2.4.0
singularity pull docker://ghcr.io/zellerlab/cayman:latest
# do something

python extract_viral_archeal_eukaryote_reads_v3.py /proj/naiss2023-23-618/FIBER_INT_METAGENOMIC_DATA/02_RAW_READS/CLEAN_READS_V1/TEST_KRAKEN_PROCESSING/M159_DC
