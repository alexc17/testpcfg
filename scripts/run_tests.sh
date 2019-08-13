#!/bin/bash
x=$1
y=$2
for (( c=$x; c<=$y; c++ ))
do
	echo "Starting $c"
   python ../testpcfg/oracle_test.py --seed 1 ../data/test/grammar${c}.pcfg ../data/test/grammar${c}.wcfg > ../data/test/grammar${c}.log &
done
wait
echo "Done"
