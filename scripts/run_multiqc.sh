#!/bin/bash

source ./config.cfg

conda activate /work/geisingerlab/conda_env/multiQC

multiqc $BASE_DIR --ignore
