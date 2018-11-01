To complile:

source set_env.sh

mkdir build

cd build

cmake ..

make

To run the program for high-dimensional mediation analysis:

./mcmc_out n c1 c2 q burnIn iterations 1 1 A_file M_file Y_file C1_file C2_file \
            beta_m.ini_file alpha_a.ini_file pi_m.ini_file pi_a.ini_file
            
            
