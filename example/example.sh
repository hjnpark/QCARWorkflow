#!/bin/bash
#SBATCH -p high
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -c 2
#SBATCH -J qcai_test
#SBATCH -t 5-00:00:00
#SBATCH --no-requeue
#--mem=16000


qcfractal-manager --fractal-uri=https://localhost:7777/  --verify False -u User1 -p 1234

#Make sure to change 'localhost' to the name of the machine that contains qcfractal-server (where you ran 'qcfractal-server start'). 
