#!/usr/bin/env python
# coding: utf-8

# In[6]:


from binIO import readin_bin, writeout_bin
#import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def get_reference(path=None, nr=None, ish=None, ich=None):
  if not path:
    path = '/home/wyw/2021_workspace/First_arrival_picking/model_field_labels'
  y_all_ = np.zeros((nr, 1), dtype=int)
  #     for ish in shot_list: # shot 
  tlabel_path = os.path.join(path + '/' + str(ich), str(ish) + '_fb' )
  #print(tlabel_path)
  if not os.path.exists(tlabel_path):
    return np.array(y_all_)
  with open(tlabel_path) as f:
    lines = f.readlines()    
    #print(f"lines, {lines}")
  for alin in lines:
    readin_c = int(float(alin.split(' ')[-1][:-1])) - 1
    readin_f = float(alin.split(' ')[0][:])
    readin_fi = int(readin_f / 6 * 3001)
    #print(f"labels {readin_c}, {readin_f}")
    y_all_[readin_c, 0] = readin_fi
  return np.array(y_all_)

# In[9]:


## cut first arrivals for S-waves
numid = 275 #sys.argv[1] # 76; 275
targtmodel = sys.argv[1] # k3 k5 k5bn k5c5 k5cat 
print("the number of shots in real data is ", numid)
print("the target network is ", targtmodel)

wkpath = '/home/wyw/2021_workspace/EXPERIENCE/new_train/model_23_gradient_clip/'
dfname = wkpath + 'real_' + str(numid) + '_' + str(targtmodel) 
os.makedirs(dfname + '_cut', exist_ok=True)
print("Done making directory: ", dfname+'_cut')

pick_path = '/home/wyw/2021_workspace/First_arrival_picking/model_field_labels'
    
for indx in range( int(numid) ):
#for indx in range( 89, 90 ):
    print(f"------ Processing Num {indx} ------")
    # path
    path_pre_s = dfname + '/' + str(indx) + '_S' 
    path_pre_p = dfname + '/' + str(indx) + '_P' 
    #path_vz = "/home/wyw/2021_workspace/B-rename_for_input/real_data_" + str(numid) + "/vz/vz_" + str(indx) + ".bin"
    path_vz = "/home/wyw/2021_workspace/B-rename_for_input/real_data/vz/vz_" + str(indx) + ".bin"
    
    # readin
    pre_s = readin_bin(path_pre_s, 0, 128, 2048)
    pre_p = readin_bin(path_pre_p, 0, 128, 2048)
    vz = readin_bin(path_vz, 0, 2048, 126)

    
    outpath_s = dfname + "_cut" + "/" + str(indx) + "_S"
    outpath_p = dfname + "_cut" + "/" + str(indx) + "_P"
    
    # cut point
#    firstpoint = np.zeros( (1, 126) , dtype=int )
    #path = '/home/wyw/2021_workspace/First_arrival_picking/model_field_labels'
    ctr = 100
    firstpoint = get_reference(path=pick_path, nr=126, ish=indx, ich='vz')
    if firstpoint[0, 0] == 0 and firstpoint[-1,0] == 0:
      print(f"no pick shot {indx}")
      for ii in range(126):
          for jj in range(2048 - ctr):
              if abs(vz[jj,ii]) > 1e-14:
                  firstpoint[ii, 0] = int(jj)
      #            print(jj, ii)
                  break
     
    # cut
#    ctr = 200
    for ii in range(128):
        #print(f"debug: first arrival point is {firstpoint[ii,0]}")
        if ii >= 126:
            pre_s[ii, 0:firstpoint[125,0] + ctr ] = 0
            pre_p[ii, 0:firstpoint[125,0] + ctr ] = 0
        else:
            pre_s[ii, 0:firstpoint[ii,0] + ctr] = 0
            pre_p[ii, 0:firstpoint[ii,0] + ctr] = 0

    #for ii in range(126):
    #    vz[0 : firstpoint[ii] + ctr, ii] = 0
    #    vx[0 : firstpoint[ii] + ctr, ii] = 0
    #    div[0 : firstpoint[ii] + ctr, ii] = 0
    #    curl[0 : firstpoint[ii] + ctr, ii] = 0
        
#     plt.imshow(pre_s.transpose(), extent=[0,5,10,0], vmax=1e-1, vmin=-1e-1)
    writeout_bin(outpath_s, pre_s, 0, 128, 2048)
    writeout_bin(outpath_p, pre_p, 0, 128, 2048)
    #writeout_bin(outpath_vz, vz, 0, 2048, 126)
    #writeout_bin(outpath_vx, vx, 0, 2048, 126)
    #writeout_bin(outpath_div, div, 0, 2048, 126)
    #writeout_bin(outpath_curl, curl, 0, 2048, 126)
    



