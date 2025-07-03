#! /bin/bash

if [ ! -e store_ramfile ]; then
    mkdir store_ramfile
fi

if [ ! -e store_epsilon ]; then
    mkdir store_epsilon
fi

for (( energy=0;energy<=300;energy++ ))
do 
  file=$(printf "%0.2f\n" $(echo "scale=2;$energy/100" | bc))
  if [ ! -e store_ramfile/RAMFILE_$file ]; then
    echo $file | genRAram610_dynamic
  
    mv RAMFILE_* store_ramfile/.
    mv EPSILON_* store_epsilon/.
  fi
done



