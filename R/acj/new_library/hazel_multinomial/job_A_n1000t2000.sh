#!/bin/bash
#BSUB -J cj_A2000[1-20]
#BSUB -n 64
#BSUB -W 2:00
#BSUB -o logs/out.%J.%I
#BSUB -e logs/err.%J.%I

module load conda
source ~/.bashrc
conda activate /usr/local/usrapps/cads/cdondim/r_catfda

cd "$LS_SUBCWD"
Rscript run_hazel_A_n1000t2000.R
