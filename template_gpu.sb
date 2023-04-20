#!/bin/bash

#SBATCH -p gpu
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --gres=gpu:nvidia_geforce_gtx_1080_ti:1
#SBATCH --mem 23000
#SBATCH --time 2-0:0:0
#SBATCH -J 3/5_name
#SBATCH -o train_name.out
#SBATCH -e train_name.err

. /home/drk/anaconda3/etc/profile.d/conda.sh
conda activate tf210

basenji_train.py
