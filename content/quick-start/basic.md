---
title: "Basic IO"
date: 2023-04-04T18:11:37+08:00
draft: true
---

## 1. read and write binary 
We use [`struct`](https://docs.python.org/3/library/struct.html)(a python package) to read binary data.

```python
#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os
import struct


def readin_bin(path=None, seek_num=None, nt=None, nr=None):
    FA = open(path, "rb")
    FA.seek(seek_num) # skip N bit to start reading in data
    out_data = np.empty((nt,nr))
    for rr in range(nr):
        for tt in range(nt):
            data = FA.read(4) # read 4 bit float
            data_float = struct.unpack("f", data)[0]
            out_data[tt][rr] = data_float
    return out_data


def writeout_bin(path=None, in_data=None, seek_num=None, nt=None, nr=None):
    FA = open(path, "wb")
    for rr in range(nr):
        for tt in range(nt):
            data = in_data[tt][rr] #FA.read(4)
            data_float = struct.pack("f", data)
            FA.write(data_float)
    FA.close()


```
## 2. convert binary to .npy | .npy to binary

```python
def npy2bin(inname=None, ouname=None, seek_num=0, nt=None, nr=None):
  indata = np.load(inname)
  writeout_bin(path=ouname, in_data=indata, seek_num=seek_num, nt=nt, nr=nr)
  return indata

def bin2npy(inname=None, ouname=None, seek_num=0, nt=None, nr=None):
  indata = readin_bin(path=inname, seek_num=seek_num, nt=nt, nr=nr)
  np.save(ouname, indata)
  return indata
```


## 3. parallel binary to .npy

find more scripts in our [github](https://github.com/weiyw16/weiyw16.github.io/blob/main/scripts/workflow/3-rename_data/template_rename/bin2npy_convert.py)