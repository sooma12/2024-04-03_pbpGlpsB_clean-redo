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

# Note, learned a nice trick for arrays, but I think you have to define the path explicitly:
# This SBATCH header would start array jobs for each file *.txt
# --array=0-$(ls *.txt | wc -l) script.sh

echo "Loading tools"
module load samtools/1.19.2

source ./config_bowtie.cfg

echo "Looking for .sam files and outputting sorted bams in $MAPPED_DIR"

# Get array of sam files
# shellcheck disable=SC2207
sams_array=($(ls -d ${MAPPED_DIR}/*.sam))

# Get specific file for this array task
current_file=${sams_array[$SLURM_ARRAY_TASK_ID-1]}

current_name=$(basename "$current_file")
current_name_no_ext="${current_name%.*}"

samtools view -bS "${current_file}" > ${MAPPED_DIR}/${current_name_no_ext}.bam
samtools sort ${MAPPED_DIR}/${current_name_no_ext}.bam -o "${MAPPED_DIR}/${current_name_no_ext}"_sorted.bam

mkdir -p $MAPPED_DIR/intermediate_files
mv ${current_file} intermediate_files/
mv "${current_name_no_ext}".bam intermediate_files/
