#!/bin/bash
#BSUB -J cj_A300
#BSUB -n 32
#BSUB -W 6:00
#BSUB -o logs/out.%J
#BSUB -e logs/err.%J

module load conda
source ~/.bashrc
conda activate /usr/local/usrapps/cads/cdondim/r_catfda

Rscript run_hazel_A_n1000t300.R
