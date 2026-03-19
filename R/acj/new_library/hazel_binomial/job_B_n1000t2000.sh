#!/bin/bash
#BSUB -J cj_bin_B2000
#BSUB -n 64
#BSUB -R "span[hosts=1]"
#BSUB -W 2:00
#BSUB -o logs/out.%J
#BSUB -e logs/err.%J

module load conda
source ~/.bashrc
conda activate /usr/local/usrapps/cads/cdondim/r_catfda

Rscript run_hazel_B_n1000t2000.R
