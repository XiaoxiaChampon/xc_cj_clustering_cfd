#!/bin/bash
#BSUB -J cj_bin_A_n100
#BSUB -n 64
#BSUB -W 2:00
#BSUB -o logs/out.%J
#BSUB -e logs/err.%J

module load conda
source ~/.bashrc
conda activate /usr/local/usrapps/cads/cdondim/r_catfda

Rscript run_hazel_A_n100t2000.R
