#!/bin/bash

DATAPATH=`pwd`/

echo in="./model_vp_$1.bin" n1=467 n2=801 data_format=naive_float > vp_raw.rsf

sfput < vp_raw.rsf \
  d1=10 d2=10 o1=0 o2=0 \
  label1=Depth label2=Distance \
  unit1=m unit2=m > vp.rsf

sfcp < vp.rsf > vp_$1.rsf

sfcp < vp.rsf > vs_raw.rsf
sfadd < vs_raw.rsf scale=0.6 mode=m > vs.rsf

sfcp < vp.rsf > rho_raw.rsf
sfadd < rho_raw.rsf scale=0 mode=m > rho_0.rsf
sfadd < rho_0.rsf add=2400 mode=a > rho.rsf

sfcat < vp.rsf vs.rsf rho.rsf > combimodel_ori.rsf 
sftransp < combimodel_ori.rsf plane=13 > newtestmodel.rsf
