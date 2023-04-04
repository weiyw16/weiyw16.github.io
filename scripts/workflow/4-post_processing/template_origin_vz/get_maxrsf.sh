#!/bin/bash

mnum=$1
#datapath=/data/wyw/make_data/case_large_train/expr-pow-2048-nusgpu-20200304/output
#datapath=/data/wyw/make_data/case_large_train/process_shotmore_model_7_pow/output
#datapath=/data/wyw/make_data/case_large_train/modeltest_model${mnum}_141/output
datapath=/home/wyw/2021_workspace/A-raw_forward_seismic/syn_test_model_${mnum}/output

for (( inn=0; inn< 140; inn++ ))
do
  ii=$(( ${inn} * 5 + 50 ))
  dataname=${datapath}/model_${mnum}_shot_${ii}_vz.bin
  echo "in=${dataname} n1=126 n2=12004 data_format=native_float" > model_${mnum}_shot_${ii}_vz_raw.rsf
  sfmath < model_${mnum}_shot_${ii}_vz_raw.rsf output="abs(input)" |\
    sfmax axis=2 max=y > model_${mnum}_shot_${ii}_vz_max.rsf
done

