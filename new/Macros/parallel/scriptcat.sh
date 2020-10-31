#!/bin/bash
numcores=12
i=1 # parallel jobs start from 1
cat mps_resultsmatrix_fake_STATSYS_after_JID1.txt > mps_resultsmatrix_fake_STATSYS_after_JID0.txt
((i = i + 1))
while [[ $i -lt $numcores ]]
do
    echo "cat $i"
    cat "mps_resultsmatrix_fake_STATSYS_after_JID$i.txt" >> mps_resultsmatrix_fake_STATSYS_after_JID0.txt
    ((i = i + 1))
done