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

# Get array of sam files
# shellcheck disable=SC2207
sams_array=($(ls -d ${DATA_DIR}/mapped/*.sam))

# Get specific file for this array task
current_file=${sams_array[$SLURM_ARRAY_TASK_ID]}

current_name=$(basename "$current_file")
current_name_no_ext="${current_name%.*}"

samtools view -bS "${current_file}" > "${current_name_no_ext}".bam
samtools sort "${current_name_no_ext}".bam -o "${current_name_no_ext}"_sorted.bam
