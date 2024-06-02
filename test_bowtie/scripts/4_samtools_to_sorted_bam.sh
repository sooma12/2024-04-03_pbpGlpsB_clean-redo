#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=samtools_to_sorted_bam
#SBATCH --time=04:00:00
#SBATCH --array=1-9%10
#SBATCH --ntasks=9
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/2024-04-03_pbpGlpsB_clean-redo/test_bowtie/logs/%x-%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/2024-04-03_pbpGlpsB_clean-redo/test_bowtie/logs/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=soo.m@northeastern.edu

echo "Loading tools"
module load samtools/1.19.2

source ./config_bowtie.cfg

for file in "${DATA_DIR}"/mapped/*.sam; do
  samtools view -bS "${file}" > "${file}".bam
  samtools sort "${file}".bam -o "${file}"_sorted.bam
done
