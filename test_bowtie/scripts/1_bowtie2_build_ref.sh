#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=bowtie2_build
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/2024-04-03_pbpGlpsB_clean-redo/test_bowtie/logs/%x-%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/2024-04-03_pbpGlpsB_clean-redo/test_bowtie/logs/%x-%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=soo.m@northeastern.edu

## Usage: sbatch 1_bowtie2_build_ref.sh

echo "loading tools for STAR genomeGenerate"
module load bowtie/2.5.2

# Load config file
source ./config_bowtie.cfg
# Relevant variables: GENOME_REF_DIR_OUT, FASTA_IN, GTF_IN
echo "Genome fasta files: $REF_CHR_FA, $REF_PAB1_FA, $REF_PAB2_FA, $REF_PAB3_FA"
echo "Basename of output: $BT2_OUT_BASE"
echo "Directory for Bowtie2 genome index: $GENOME_INDEX_DIR"

bowtie2-build -f --threads 4 $REF_CHR_FA,$REF_PAB1_FA,$REF_PAB2_FA,$REF_PAB3_FA $BT2_OUT_BASE

# mkdir -p $GENOME_INDEX_DIR
# mv *.bt2 $GENOME_INDEX_DIR
