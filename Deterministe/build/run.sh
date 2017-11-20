#! /bin/sh

expe_type=2
nb_segs=100
mu=1
nb_pts_mu=10
sigma_at=1
sigma_a_arg1=1
sigma_st=1
sigma_s_arg1=2
sourcet=1
source_arg1=1
epsilon=1e-3
output_style=3
fname=output

echo $expe_type $nb_segs $mu $nb_pts_mu $sigma_at $sigma_a_arg1 $sigma_st $sigma_s_arg1 $epsilon $sourcet $source_arg1 $output_style $fname 

make && ./solver_deter $expe_type $nb_segs $mu $nb_pts_mu $sigma_at $sigma_a_arg1 $sigma_st $sigma_s_arg1 $epsilon $sourcet $source_arg1 $output_style $fname
