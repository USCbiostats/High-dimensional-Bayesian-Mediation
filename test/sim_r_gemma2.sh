#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --job-name=PMY
#SBATCH --mem=10000
#SBATCH --array=1-1000%250
#SBATCH --workdir=/net/wonderland/home/yanys/Bayesian_mediation/RealSimulation/BF_simulation/permutation 

bash

let k=0

for((m=1;m<=1;m++)); do
  for((rep=1;rep<=1000;rep++)); do
    for((i=1;i<=1;i++)); do
    
      let k=${k}+1
      if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
      let start=(${i}-1)*10+1
      let end=${i}*10
      for((permute=1;permute<=100;permute++)); do
        echo ${permute}
        if [ -f results_pM_1${rep}_${permute}.txt ]; then
          echo 'find files';
        fi
        if [ ! -f results_pM_1${rep}_${permute}.txt ]; then
        MCMC=/net/wonderland/home/yanys/Bayesian_mediation/RealSimulation/BF_simulation/Bayesian-Mediation-Analysis/build/bin/mcmc_M.x
        Rscript permuteM_rep1.R ${m} ${rep} ${permute}
        ${MCMC} 1000 1 1 100 50000 100000 ${m} ${rep} ${permute} Y_${m}${rep}.txt permuteM_${m}${rep}_${permute}.txt A_${m}${rep}.txt C1.txt C1.txt beta_m_${m}${rep}_${permute}.txt alpha_a_${m}${rep}_${permute}.txt pi_m_${m}${rep}_${permute}.txt pi_a_${m}${rep}_${permute}.txt  
        fi
        
        if [ -f results_pY_1${rep}_${permute}.txt ]; then
          echo 'find files';
        fi
        if [ ! -f results_pY_1${rep}_${permute}.txt ]; then
        MCMC=/net/wonderland/home/yanys/Bayesian_mediation/RealSimulation/BF_simulation/Bayesian-Mediation-Analysis/build/bin/mcmc.x
        Rscript permuteY_rep1.R ${m} ${rep} ${permute}
        ${MCMC} 1000 1 1 100 50000 100000 ${m} ${rep} ${permute} permuteY_${m}${rep}_${permute}.txt M_${m}${rep}.txt A_${m}${rep}.txt C1.txt C1.txt beta_m_${m}${rep}_${permute}.txt alpha_a_${m}${rep}_${permute}.txt pi_m_${m}${rep}_${permute}.txt pi_a_${m}${rep}_${permute}.txt  
        fi
        
        rm -f permuteM_${m}${rep}_${permute}.txt
        rm -f permuteY_${m}${rep}_${permute}.txt
        rm -f pi_m_${m}${rep}_${permute}.txt
        rm -f pi_a_${m}${rep}_${permute}.txt
        rm -f beta_m_${m}${rep}_${permute}.txt
        rm -f alpha_a_${m}${rep}_${permute}.txt
        rm -f Y_${m}${rep}.txt
        rm -f M_${m}${rep}.txt
        rm -f A_${m}${rep}.txt
      
      done
      fi
    done
  done
done

