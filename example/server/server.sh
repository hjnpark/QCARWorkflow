#!/bin/bash
#SBATCH -p high
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 1
#SBATCH -J qcai_server
#SBATCH -t 10-00:00:00
#SBATCH --no-requeue
#--mem=16000


qcfractal-server start --base-folder ~/.qca/qcfractal 

