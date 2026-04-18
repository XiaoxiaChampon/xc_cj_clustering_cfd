#!/bin/bash
#BSUB -J cj_bin_A_n100
#BSUB -n 64
#BSUB -W 2:00
#BSUB -o logs/out.%J
#BSUB -e logs/err.%J

module load conda
source ~/.bashrc
conda activate /path/to/your/r_catfda_env  # adjust to your HPC environment

Rscript run_hazel_A_n100t2000.R
