## convert bin to npy
## made by Yanwen Wei @ 2021-4-26
##
import os
import numpy as np
import struct
from joblib import Parallel, delayed
import sys
NT = 2048 #int(sys.argv[1])
NR = 126 #int(sys.argv[2])

## function
def readin_bin(path, seek_num=0, nt=None, nr=None):
        FA = open(path, "rb")
        FA.seek(seek_num)
        out_data = np.empty((nt, nr))
        for rr in range(nr):
            for tt in range(nt):
                data = FA.read(4)
                data_float = struct.unpack("f", data)[0]
                out_data[tt][rr] = data_float
        return out_data
      
def sub_save(this_index=None, nt=None, nr=None, mod=None):
  this_path_mix = f"./{mod}/{mod}_" + str(this_index) + '.bin'
  that_path_mix = f"./{mod}/{mod}_" + str(this_index) + '.npy'
  tmp_in = readin_bin(this_path_mix, 0, nt, nr)
  print(f"{this_index}, {(np.abs(tmp_in)).max()}")
  np.save( that_path_mix, tmp_in)
  return 


def save_data(this_index=None):
  sub_save(this_index=this_index, nt=NT, nr=NR, mod='vz')
  sub_save(this_index=this_index, nt=NT, nr=NR, mod='vx')
  sub_save(this_index=this_index, nt=NT, nr=NR, mod='div')
  sub_save(this_index=this_index, nt=NT, nr=NR, mod='curl')

  return 

def main():

  all_mix = []
  od = os.getcwd()
  path = od + '/vz'      # 输入文件夹地址
  files = os.listdir(path)   # 读入文件夹
  for ii in files:
    if '.bin' in ii:
      all_mix.append(ii)
      
  num_png = len(all_mix)
  print(f"There are {num_png} files.")

  results = (
      Parallel(n_jobs=20)
      (delayed(save_data)(ii)
        for ii in range(num_png))
      )

if __name__ == '__main__':
  main()

