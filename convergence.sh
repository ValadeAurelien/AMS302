#! /bin/sh

N_exp_min=2
N_exp_max=7
N_pts=14
stype=2 #type source (1->cste, 2->delta(0))
file=err_N_$stype

echo > $file

for step in $(seq 1 $N_pts)
do
    echo $step
    exp=$(calc "$step*($N_exp_max-$N_exp_min)/$N_pts+$N_exp_min" | sed 's/~//g')
    N=$(calc "floor(10**($exp))" | sed 's/~//g')
    ./TP1 1 1 $N 100 $stype 3 tmp
    err=$(awk '/diff/ {print $4}' tmp)
    echo $exp $N $err
    echo $N $err >> $file
done
