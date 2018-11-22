#!/bin/bash

#
# Generate test data for special functions
#

set -e

mkdir -p data/special
g++ -o gen_special_test gen_special_test.cpp
N=1000
for fun in betaInt beta ibetaInt ibeta choose erf erfc gamma igamma;
do
    echo "generating data for special/${fun}"
    r=$(./mkSeed ${fun})
    ./gen_special_test ${fun} ${N} ${r} > data/special/${fun}.yaml
done
