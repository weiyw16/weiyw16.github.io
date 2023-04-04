#!/usr/bin/env python
# coding: utf-8

# In[6]:


from binIO import readin_bin, writeout_bin
#import matplotlib.pyplot as plt
import numpy as np
import os
import sys

########################################################
orginmodel = sys.argv[1] 
targtmodel = orginmodel #sys.argv[2] 
print("the original model num is ", orginmodel)
print("the target model num is ", targtmodel)

os.makedirs('output_rmp', exist_ok=True)
inpre = 'output/'
oupre = 'output_rmp/'
pick_path = '/home/wyw/2021_workspace/First_arrival_picking/ray_trace_calculate/tomo/proj_8/caltime.tt'
nt = 12004
ns = 71
########################################################


def get_reference(path=None, nr=None, ish=None, ich=None):
  if not path:
    #path = '/home/wyw/2021_workspace/First_arrival_picking/model_field_labels'
    path = '/home/wyw/2021_workspace/First_arrival_picking/ray_trace_calculate/tomo/proj_8/caltime.tt'
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
    readin_fi = int(readin_f / 4 * 12004)
    #print(f"labels {readin_c}, {readin_f}")
    y_all_[readin_c, 0] = readin_fi
  return np.array(y_all_)


def get_calculation_tomo(path=None):
  # read tt and plot
  rsid=[]; rnr=[]; rsx=[]; rsz=[]
  Rrid=[]; Rrx=[]; Rrz=[]; Rtt=[]
  #     f = open('./caltime.tt','r')
  f = open(path,'r')
  tthead = f.readline()
  rns = int(tthead.split()[0])
  for si in range(0,rns):
    shead1 = f.readline()
    rsid.append(int(shead1.split()[0]))
    rnr.append(int(shead1.split()[1]))
    shead2 = f.readline()
    rsx.append(float(shead2.split()[0]))
    rsz.append(float(shead2.split()[1]))
    rrid=[]; rrx=[]; rrz = []; rtt=[]
    for ri in range(0,rnr[si]):
      ttime = f.readline()
      rrid.append(int(ttime.split()[0]))
      rrx.append(float(ttime.split()[1]))
      rrz.append(float(ttime.split()[2]))
      rrt.append(float(ttime.split()[3]))
    Rrid.append(rrid)
    Rrx.append(rrx)
    Rrz.append(rrz)
    Rtt.append(rtt)
  return Rtt


for indx in range(ns):
    # path
    path_ref_s = inpre + "model_"     + str(orginmodel) + "_shot_" + str(indx) + "_curl.bin"
    path_ref_p = inpre + "model_"     + str(orginmodel) + "_shot_" + str(indx) + "_div.bin"
    path_vz    = inpre + "model_"     + str(orginmodel) + "_shot_" + str(indx) + "_vz.bin"
    path_vx    = inpre + "model_"     + str(orginmodel) + "_shot_" + str(indx) + "_vx.bin"
    
    # readin
    #pre_s = readin_bin(path_pre_s, 0, 128, nt)
    #pre_p = readin_bin(path_pre_p, 0, 128, nt)
    div = readin_bin(path_ref_p, 0, nt, 126)
    curl = readin_bin(path_ref_s, 0, nt, 126)
    vz = readin_bin(path_vz, 0, nt, 126)
    vx = readin_bin(path_vx, 0, nt, 126)

    
    #outpath_s = "/home/wyw/workspace/EXPERIMENTS/job/job_bgp_2to1/newbin_2to2/testmodel_"     + str(targtmodel) + "_141/" + str(indx) + "_S"
    #outpath_p = "/home/wyw/workspace/EXPERIMENTS/job/job_bgp_2to1/newbin_2to2/testmodel_"     + str(targtmodel) + "_141/" + str(indx) + "_P"
    outpath_div  = oupre + "model_"     + str(targtmodel) + "_shot_" + str(indx) + "_div.bin"
    outpath_curl = oupre + "model_"     + str(targtmodel) + "_shot_" + str(indx) + "_curl.bin"
    outpath_vz   = oupre + "model_"     + str(targtmodel) + "_shot_" + str(indx) + "_vz.bin"
    outpath_vx   = oupre + "model_"     + str(targtmodel) + "_shot_" + str(indx) + "_vx.bin"

    ctr = 100
    #firstpoint = get_reference(path=pick_path, nr=126, ish=indx, ich='vz')
    firstpoint = get_calculation_tomo(path=pick_path)[indx]
    print(firstpoint.shape)
    if firstpoint[0, 0] == 0 and firstpoint[-1,0] == 0:
      print(f"no pick shot {indx}")
      for ii in range(126):
          for jj in range(nt - ctr):
              if abs(vz[jj,ii]) > 1e-14:
                  firstpoint[ii, 0] = int(jj)
      #            print(jj, ii)
                  break
    
    # cut
    for ii in range(126):
        vz[0 : firstpoint[ii] + ctr, ii] = 0
        vx[0 : firstpoint[ii] + ctr, ii] = 0
        div[0 : firstpoint[ii] + ctr, ii] = 0
        curl[0 : firstpoint[ii] + ctr, ii] = 0
        
#     plt.imshow(pre_s.transpose(), extent=[0,5,10,0], vmax=1e-1, vmin=-1e-1)
    #writeout_bin(outpath_s, pre_s, 0, 128, nt)
    #writeout_bin(outpath_p, pre_p, 0, 128, nt)
    writeout_bin(outpath_vz, vz, 0, nt, 126)
    writeout_bin(outpath_vx, vx, 0, nt, 126)
    writeout_bin(outpath_div, div, 0, nt, 126)
    writeout_bin(outpath_curl, curl, 0, nt, 126)
    



