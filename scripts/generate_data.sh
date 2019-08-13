#!/bin/bash
seed=1
n=2
sampler="../../syntheticpcfg/syntheticpcfg/sample_grammar.py"
for br in 20 30 40 50 60 70 80
do
    direct="../data/test$br"
    mkdir $direct
    python $sampler --seed $seed --numbergrammars $n --binaryproductions $br "${direct}/grammar%d.pcfg"
done
echo "Done"
