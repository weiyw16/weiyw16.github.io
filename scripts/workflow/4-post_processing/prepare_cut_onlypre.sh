#!/bin/bash 

orinum=$1
tarnum=$2

cp -Lr template_p ./model${tarnum}_p
#cp -Lr template_origin_p ./model${tarnum}_origin_p
cp -Lr template_s ./model${tarnum}_s
#cp -Lr template_origin_s ./model${tarnum}_origin_s

cd ./model${tarnum}_p
bash ./get_maxrsf.sh ${orinum}
bash ./restore_pre.sh ${orinum} ${tarnum}
cd ..
#cd ./model${tarnum}_origin_p
#bash ./get_maxrsf.sh ${orinum}
#bash ./restore_pre.sh ${orinum} ${tarnum}
#cd ..
cd ./model${tarnum}_s
bash ./get_maxrsf.sh ${orinum}
bash ./restore_pre.sh ${orinum} ${tarnum}
cd ..
#cd ./model${tarnum}_origin_s
#bash ./get_maxrsf.sh ${orinum}
#bash ./restore_pre.sh ${orinum} ${tarnum}
#cd ..

