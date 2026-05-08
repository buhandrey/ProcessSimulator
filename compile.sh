#!/bin/bash

for file in ./*.cpp; do
    fullname=${file##*/};
    basename=${fullname%.cpp};
    numcol=$(awk '{print NF}' ${file} | sort -nu | tail -n 1);
    echo "${file} (${basename})";
    g++ -o ${basename} ${fullname} -O3 -fopenmp;
done

exit
