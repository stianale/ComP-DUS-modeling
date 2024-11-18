#!/bin/bash

for num in {1..10};
do
newdir=$(echo $1 | sed "s/ComP_scripts\///g" | sed "s/.py//g")
mkdir ComP_outputs/"${newdir}"_"${num}"
python $1
cp /tmp/outputs/* ComP_outputs/"${newdir}"_"${num}"
done
