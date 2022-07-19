#!/bin/bash

filedr=/projects/dame8201/datasets/111_test/single-cell-models/


##for i in {1..1000}; do  if [ -f home_cell_${i}.mat ]; then echo ""; else echo ${i} >> missing_flz.txt; fi; done;


for i in {1..1000}
do 
if ! [ -f home_cell_${i}.mat ]
then
echo "file ${i} not found"
fi
done
