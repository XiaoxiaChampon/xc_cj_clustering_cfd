#!/bin/bash
#BSUB -J cj_A_n100[1-10]
#BSUB -n 32
#BSUB -W 6:00
#BSUB -o logs/out.%J.%I
#BSUB -e logs/err.%J.%I

module load conda
source ~/.bashrc
conda activate /usr/local/usrapps/cads/cdondim/r_catfda

Rscript run_hazel_A_n100t2000.R
