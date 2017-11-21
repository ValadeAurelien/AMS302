#! /bin/sh

##############
### HELPER ###
#########################################################################
# expe_type     = 1 (pas de diff) , 2 (diff), 3 (DSA)                   #
# nb_segs       = partition de l'espace				        #
# mu            = valeur de mu (quand utilisée, sinon inutile)	        #
# nb_pts_mu     = partition des mus (quand utilisées, sinon inutile)	#
# sigma_at      = 1 (constante), 3 (marches)				#
# sigma_a_arg1  = valeur de sigma absorbtion				#
# sigma_st      = 1 (constante), 3 (marches)				#
# sigma_s_arg1  = valeur de sigma diffusion				#
# sourcet       = 1 (constante), 2 (delta(0))			        #
# source_arg1   = valeur de la source				        #
# epsilon       = critère d'arret de la boucle quand diffusion	        #
# output_style  = 1 (plot), 2 (file), 3 (les deux), autre (rien)	#
# fname         = nom du fichier en sortie 				#
#########################################################################


########################
### CHOIX EXPERIENCE ###
########################
EXPE=bac_a_sable
#EXPE=two_steps
#EXPE=scattering
#EXPE=loop_nb_pts_mu
#EXPE=epsilon_abs
#EXPE=loop_nb_segs

function do_it {
    echo 
    echo "EXPE $1
         expe_type $2
         nb_segs $3
         mu $4
         nb_pts_mu $5
         sigma_at $6
         sigma_a_arg1 $7
         sigma_st $8
         sigma_s_arg1 $9
         epsilon ${10}
         sourcet ${11}
         source_arg1 ${12}
         output_style ${13}
         fname ${14}" | column -t -c 4 | sed 's/ /./g'
    echo $1 : $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} 
    make && ./solver_deter $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} 
}

if [[ $EXPE == "bac_a_sable" ]]
then
    expe_type=3
    nb_segs=1000
    mu=1
    nb_pts_mu=100
    sigma_at=3
    sigma_a_arg1=1
    sigma_st=1
    sigma_s_arg1=1000
    sourcet=1
    source_arg1=1
    epsilon=1e-3
    output_style=1
    fname=./datas/output_${expe_type}_TS_${sigma_a_arg1}_${sourcet}
    
    do_it $EXPE $expe_type $nb_segs $mu $nb_pts_mu $sigma_at $sigma_a_arg1 $sigma_st $sigma_s_arg1 $epsilon $sourcet $source_arg1 $output_style $fname 
elif [[ $EXPE == "two_steps" ]]
then
    expe_type=1
    nb_segs=1000
    mu=1
    nb_pts_mu=10
    sigma_at=3
    sigma_a_arg1=3
    sigma_st=1
    sigma_s_arg1=1
    sourcet=2
    source_arg1=1
    epsilon=1e-3
    output_style=3
    fname=./datas/output_two_steps_${expe_type}_${mu}_${sigma_a_arg1}_${sourcet}
    
    do_it $EXPE $expe_type $nb_segs $mu $nb_pts_mu $sigma_at $sigma_a_arg1 $sigma_st $sigma_s_arg1 $epsilon $sourcet $source_arg1 $output_style $fname
elif [[ $EXPE == "scattering" ]]
then
    expe_type=2
    nb_segs=1000
    mu=1  # pas d'importance ici
    nb_pts_mu=100
    sigma_at=2
    sigma_a_arg1=5
    sigma_st=1
    sigma_s_arg1=1
    sourcet=2
    source_arg1=1
    epsilon=1e-3
    output_style=3
    fname=./datas/output_scattering_${mu}_${sigma_a_arg1}_${sourcet}
    
    do_it $EXPE $expe_type $nb_segs $mu $nb_pts_mu $sigma_at $sigma_a_arg1 $sigma_st $sigma_s_arg1 $epsilon $sourcet $source_arg1 $output_style $fname
elif [[ $EXPE == "loop_nb_pts_mu" ]]
then
    expe_type=2
    nb_segs=1000
    mu=1  # pas d'importance ici
    sigma_at=1
    sigma_a_arg1=2
    sigma_st=1
    sigma_s_arg1=1
    sourcet=1
    source_arg1=1
    epsilon=1e-3
    output_style=3

    dirname=./datas/loop_nb_pts_mu
    mkdir -p $dirname 2>/dev/null
    
    for nb_pts_mu in $(seq 2 2 16)
    do
	echo $nb_pts_mu
	fname=$dirname/output_${sourcet}_${nb_pts_mu}
	do_it $EXPE $expe_type $nb_segs $mu $nb_pts_mu $sigma_at $sigma_a_arg1 $sigma_st $sigma_s_arg1 $epsilon $sourcet $source_arg1 $output_style $fname
    done
elif [[ $EXPE == "epsilon_abs" ]]
then

    nb_segs=200
    mu=1
    nb_pts_mu=50
    sigma_at=1
    sigma_st=1
    sourcet=1
    epsilon=1e-3
    output_style=0
    fname=output_eps

    SIGMA_A=3
    SIGMA_T=4
    fall=./datas/output_eps_nb_ite
    echo "#nb iterations for sigma_a = epsilon" > $fall
    for EPSILON in 0.1 0.01 0.001
    do
	sigma_a_arg1=$(calc "$SIGMA_A*$EPSILON")
	sigma_s_arg1=$(calc "$SIGMA_T/$EPSILON-$SIGMA_A*$EPSILON")
	source_arg1=$EPSILON
	echo >> $fall
	echo \#$EPSILON >> $fall
	
	expe_type=2
	do_it $EXPE $expe_type $nb_segs $mu $nb_pts_mu $sigma_at $sigma_a_arg1 $sigma_st $sigma_s_arg1 $epsilon $sourcet $source_arg1 $output_style $fname >> $fall

	expe_type=3
	do_it $EXPE $expe_type $nb_segs $mu $nb_pts_mu $sigma_at $sigma_a_arg1 $sigma_st $sigma_s_arg1 $epsilon $sourcet $source_arg1 $output_style $fname >> $fall
    done
elif [[ $EXPE == "loop_nb_segs" ]]
then
    expe_type=1
    mu=.5
    nb_pts_mu=50
    sigma_at=1
    sigma_a_arg1=2
    sigma_st=1
    sigma_s_arg1=1
    sourcet=2
    source_arg1=1
    epsilon=1e-3
    output_style=2

    dirname=./datas/loop_nb_segs
    mkdir -p $dirname 2>/dev/null
    fall=$dirname/output_${sourcet}_all
    
    echo "#Error as a function of the number of points" > $fall
    for nb_segs in $(seq 50 50 1000)
    do
	echo $nb_segs
	fname=$dirname/output_${sourcet}_${nb_segs}
	do_it $EXPE $expe_type $nb_segs $mu $nb_pts_mu $sigma_at $sigma_a_arg1 $sigma_st $sigma_s_arg1 $epsilon $sourcet $source_arg1 $output_style $fname
	diff=$(awk '/diff/ {print $4}' $fname)
	echo $nb_segs $diff >> $fall
    done
fi
	       
	       
	       
	       
	       
	       
	       
	       
	       
