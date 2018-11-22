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

mkdir -p data/dist
g++ -o gen_dist_test gen_dist_test.cpp
N=1000
for fun in beta binom chisq gamma geom hyper norm pois stud
do
    echo "generating data for special/${fun}"
    r=$(./mkSeed ${fun})
    ./gen_dist_test ${fun} ${N} ${r} > data/dist/${fun}.yaml
done
