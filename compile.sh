#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --job-name=permutation1
#SBATCH --mem=1000
#SBATCH --array=1-1
#SBATCH --workdir=/net/wonderland/home/yanys/Bayesian_mediation/RealSimulation/BF_simulation/Bayesian-Mediation-Analysis/

bash

lsb_release -a
source set_env.sh
mkdir build
cd build
cmake ..
make



