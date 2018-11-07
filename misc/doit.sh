#!/bin/bash

g++ -o gen_special_test gen_special_test.cpp

N=1000
for fun in beta ibeta choose erf erfc gamma igamma;
do
    r=$RANDOM
    ./gen_special_test ${fun} python ${N} ${r} > ../python/iid/test/data_${fun}.py
done
