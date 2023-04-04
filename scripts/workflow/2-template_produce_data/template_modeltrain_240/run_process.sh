#!/bin/bash
model_num=0
gpuid=5
mkdir output
exepath=/home/wyw/2021_workspace/vsp-ps-separation-workspace/1-forward_modeling_src/bin
ln -bs /home/wyw/2021_workspace/vsp-ps-separation-workspace/0-velocity_models/newtestmodel.rsf@ ./modelfile_${model_num}.bin

for((ssnum=50; ssnum < 750; ssnum+=401));
do

for (( spos=${ssnum}; spos < ${ssnum}+401; spos+=50 ));
do
  echo "processing MODEL ${model_num} SOURCE ${spos}"
  echo -e "801 467 10.0 10.0\n5 ${spos} 5 1\n15 50 16 50 2 126\n15.0 1000.0 0.0005 12004\n120 120\n1\n0 200\n0 0\n" > ParamInput_${spos}.txt
  ${exepath}/fd2d_snap_mprocess_multirec ${model_num} ${spos} ${gpuid} & 
done
wait
done
