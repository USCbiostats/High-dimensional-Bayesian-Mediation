#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --job-name=PMY
#SBATCH --mem=10000
#SBATCH --array=1-1000%250
#SBATCH --workdir=/net/wonderland/home/yanys/Bayesian_mediation/RealSimulation/BF_simulation/permutation 

bash


MCMC_M=/net/wonderland/home/yanys/Bayesian_mediation/RealSimulation/BF_simulation/Bayesian-Mediation-Analysis/build/bin/mcmc_M.x
MCMC_0=/net/wonderland/home/yanys/Bayesian_mediation/RealSimulation/BF_simulation/Bayesian-Mediation-Analysis/build/bin/mcmc.x

n=1000
w1=1
w2=1
q=100
bunrin=50000
niter=100000



let k=0

for((m=1;m<=1;m++)); do
  for((rep=1;rep<=1000;rep++)); do
    for((i=1;i<=1;i++)); do
    
      let k=${k}+1

      if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
        let start=(${i}-1)*10+1
        let end=${i}*10

        for((permute=1;permute<=100;permute++)); 
        do
          echo ${permute}
          Rscript permuteM_rep1.R ${m} ${rep} ${permute}
          Rscript permuteY_rep1.R ${m} ${rep} ${permute}
        done

        ${MCMC_M} $n $w1 $w2 $q $burnin $niter ${m} ${rep} ${permute} \
          Y_${m}${rep}                                                \
          permuteM_${m}${rep}                                         \
          A_${m}${rep}                                                \
          C1                                                          \
          C1                                                          \
          beta_m_${m}${rep}                                           \
          alpha_a_${m}${rep}                                          \
          pi_m_${m}${rep}                                             \
          pi_a_${m}${rep}   
          
        ${MCMC_0} $n $w1 $w2 $q $burnin $niter ${m} ${rep} ${permute} \
          permuteY_${m}${rep}                                         \
          M_${m}${rep}                                                \
          A_${m}${rep}                                                \
          C1                                                          \
          C1                                                          \
          beta_m_${m}${rep}                                           \
          alpha_a_${m}${rep}                                          \
          pi_m_${m}${rep}                                             \
          pi_a_${m}${rep}  
          
        for((permute=1;permute<=100;permute++)); 
        do
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

