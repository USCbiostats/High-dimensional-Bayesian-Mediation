#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --job-name=PMY
#SBATCH --mem=10000
#SBATCH --array=1-10
#SBATCH --workdir=/net/mulan/yanys/data/MESA/permutation

bash

export OMP_NUM_THREADS=24

MCMC_0=/net/wonderland/home/yanys/Bayesian_mediation/RealSimulation/BF_simulation/Bayesian-Mediation-Analysis/build/bin/mcmc_rd_omp.x
MCMC_Y=/net/wonderland/home/yanys/Bayesian_mediation/RealSimulation/BF_simulation/Bayesian-Mediation-Analysis/build/bin/mcmc_rdY_omp.x


n=1224
w1=5
w2=9
q=10000
burnin=50000
niter=100000
Y=/net/mulan/yanys/data/MESA/NGB/dbp5c_ngb_st.txt
M=/net/mulan/yanys/data/MESA/NGB/dbp5c_neighbor2_M_st.txt
A=/net/mulan/yanys/data/MESA/NGB/A_SOCENV_CEB_dbp5c_st
C1=/net/mulan/yanys/data/MESA/NGB/dbp5c_ngb_cov1_st
C2=/net/mulan/yanys/data/MESA/NGB/dbp5c_ngb_cov2_st

let k=0

for((m=1;m<=1;m++)); do
  for((rep=1;rep<=1;rep++)); do
    for((i=1;i<=10;i++)); do
    
      let k=${k}+1

      if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
        let start=(${i}-1)*10+1
        let end=${i}*10
        let permute=101

        #echo ${permute}
        #Rscript permuteM_rep1.R ${m} ${rep} ${permute} ${Y} ${M}

        #${MCMC_Y} $n $w1 $w2 $q $burnin $niter 1 ${permute}           \
        #  permuteY_${m}${rep}                                         \
        #  M_${m}${rep}                                                \
        #  ${A}                                                        \
        #  ${C1}                                                       \
        #  ${C2}                                                       \
        #  beta_m_${m}${rep}                                           \
        #  alpha_a_${m}${rep}                                          \
        #  pi_m_${m}${rep}                                             \
        #  pi_a_${m}${rep}   
        
        echo ${start}
        ${MCMC_0} $n $w1 $w2 $q $burnin $niter 2 ${start}           \
          Y_${m}${rep}                                                \
          M_${m}${rep}                                                \
          ${A}                                                        \
          ${C1}                                                       \
          ${C2}                                                       \
          beta_m_${m}${rep}                                           \
          alpha_a_${m}${rep}                                          \
          pi_m_${m}${rep}                                             \
          pi_a_${m}${rep}  
          
        #for((permute=1;permute<=100;permute++)); 
        #do
        #  rm -f permuteM_${m}${rep}_${permute}.txt
        #  rm -f permuteY_${m}${rep}_${permute}.txt
        #  rm -f pi_m_${m}${rep}_${permute}.txt
        #  rm -f pi_a_${m}${rep}_${permute}.txt
        #  rm -f beta_m_${m}${rep}_${permute}.txt
        #  rm -f alpha_a_${m}${rep}_${permute}.txt
        #  rm -f Y_${m}${rep}.txt
        #  rm -f M_${m}${rep}.txt
        #  rm -f A_${m}${rep}.txt
        #done
      fi
    done
  done
done

