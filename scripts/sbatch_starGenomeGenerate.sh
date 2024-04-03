#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=starGenomeGenerate
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/2024-04-03_pbpGlpsB_clean-redo/logs/%x-%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/2024-04-03_pbpGlpsB_clean-redo/logs/%x-%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=soo.m@northeastern.edu

## Usage: sbatch sbatch_starGenomeGenerate.sh

echo "loading tools for STAR genomeGenerate"
module load star/2.7.11a

# Load config file
source ./config.cfg
# Relevant variables: GENOME_REF_DIR_OUT, FASTA_IN, GTF_IN
echo "genome fasta file used: $FASTA_IN"
echo "genome gtf file used: $GTF_IN"
echo "genome reference (output) location: $GENOME_REF_DIR"

NTHREADS=4

# STAR requires the output directory be pre-made
mkdir -p $GENOME_REF_DIR

# STAR time
# Program recommended `--genomeSAindexNbases 9` after running with default value 14
STAR --runMode genomeGenerate \
--genomeDir $GENOME_REF_DIR \
--genomeFastaFiles $FASTA_IN \
--sjdbGTFfile $GTF_IN \
--sjdbGTFfeatureExon gene \
--genomeSAindexNbases 9 \
--runThreadN $NTHREADS
# --sjdbGTFtagExonParentTranscript Parent
