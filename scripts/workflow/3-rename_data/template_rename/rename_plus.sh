#!/bin/bash 

datapath=/home/wyw/2021_workspace/A-raw_forward_seismic/syn_test_model_$1/seismo_pow
datapath_target=`pwd`
model_range=1
shot_range=141

mkdir ${datapath_target}/vz
mkdir ${datapath_target}/vx
mkdir ${datapath_target}/div
mkdir ${datapath_target}/curl

for (( ii=0; ii<${model_range}; ii++ ))
do
  for (( jnn=0; jnn<${shot_range}; jnn++ ))
  do
    ii=$1
    jj=$(( 50 + 5 * ${jnn})) 
    #this_index=$(( ${ii}*${shot_range} + ${jj} ))
    this_index=$((  ${jnn} ))
    echo ${this_index}
    tarvz=vzmodel_${ii}_shot_${jj}_vz_afterpowmodel_${ii}_shot_${jj}_vz_plot_ready.rsf@
    tarvx=vxmodel_${ii}_shot_${jj}_vx_afterpowmodel_${ii}_shot_${jj}_vx_plot_ready.rsf@
    tardiv=divmodel_${ii}_shot_${jj}_div_afterpowmodel_${ii}_shot_${jj}_div_plot_ready.rsf@
    tarcurl=curlmodel_${ii}_shot_${jj}_curl_afterpowmodel_${ii}_shot_${jj}_curl_plot_ready.rsf@
    cp ${datapath}/${tarvz} ${datapath_target}/vz/vz_${this_index}.bin
    cp ${datapath}/${tardiv} ${datapath_target}/div/div_${this_index}.bin
    cp ${datapath}/${tarvx} ${datapath_target}/vx/vx_${this_index}.bin
    cp ${datapath}/${tarcurl} ${datapath_target}/curl/curl_${this_index}.bin
  done
done

