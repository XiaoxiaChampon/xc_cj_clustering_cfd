#!/bin/bash
#BSUB -J cj_bin_B2000
#BSUB -n 64
#BSUB -R "span[hosts=1]"
#BSUB -W 2:00
#BSUB -o logs/out.%J
#BSUB -e logs/err.%J

module load conda
source ~/.bashrc
conda activate /path/to/your/r_catfda_env  # adjust to your HPC environment

Rscript run_hazel_B_n1000t2000.R
