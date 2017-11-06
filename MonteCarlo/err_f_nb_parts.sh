#! /bin/sh

for i in {2..15}
do
    nb_parts=$(echo "sqrt(10)^$i" | bc)
    ./TP1 1 1 $nb_parts 100 1 3 tmp
    dist=$(tail -n 1 tmp | cut -d ' ' -f 4)
    echo $nb_parts $dist
done
