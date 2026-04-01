#!/bin/bash
#BSUB -J cj_bin_A750
#BSUB -n 64
#BSUB -R "span[hosts=1]"
#BSUB -W 2:00
#BSUB -o logs/out.%J.%I
#BSUB -e logs/err.%J.%I

module load conda
source ~/.bashrc
conda activate /usr/local/usrapps/cads/cdondim/r_catfda

cd "$LS_SUBCWD"
Rscript run_hazel_A_n1000t750.R
