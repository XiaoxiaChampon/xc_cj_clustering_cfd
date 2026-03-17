#!/bin/bash
#BSUB -J cj_A2000_b1
#BSUB -n 64
#BSUB -W 3:00
#BSUB -o logs/out.%J
#BSUB -e logs/err.%J

module load conda
source ~/.bashrc
conda activate /usr/local/usrapps/cads/cdondim/r_catfda

cd "$LS_SUBCWD"
export LSB_JOBINDEX=1
Rscript run_hazel_A_n1000t2000.R